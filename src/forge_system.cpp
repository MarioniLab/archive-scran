#include "scran.h"

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <class V, class M>
SEXP subset_and_divide_internal(const M in, SEXP inmat, SEXP row_subset, SEXP col_subset) {
    // Checking subset vectors
    auto rsubout=check_subset_vector(row_subset, in->get_nrow());
    const size_t rslen=rsubout.size();
    auto csubout=check_subset_vector(col_subset, in->get_ncol());
    const size_t cslen=csubout.size();

    /* Checking which rows are non-zero, and which are to be retained.
     * This is done in C++ so as to avoid needing to create the normalized expression
     * matrix and then subset it (i.e., two sets of writes).
     */
    V incoming(in->get_nrow());
    std::deque<size_t> to_retain;
    size_t start_row=0, end_row=0;
    {
        // Computing the row sums.
        V combined(rslen);
        for (const auto& c : csubout) { 
            auto inIt=in->get_const_col(c, incoming.begin());
            auto coIt=combined.begin();
            for (auto rsIt=rsubout.begin(); rsIt!=rsubout.end(); ++rsIt, ++coIt) {
                (*coIt)+=*(inIt + *rsIt);
            }
        }
        
        // Storing the indices of elements to be retained (w.r.t. the original matrix, and to the row subset).
        auto coIt = combined.begin();
        auto rsIt = rsubout.begin(); 
        for (size_t rs=0; rs<rslen; ++coIt, ++rsIt, ++rs) {
            if (*coIt >= 0.00000001) {
                to_retain.push_back(*rsIt);
            }
        }

        if (!to_retain.empty()) {
            // Cutting out extraction costs for unneeded start/end elements.
            start_row=*std::min_element(to_retain.begin(), to_retain.end());
            end_row=*std::max_element(to_retain.begin(), to_retain.end())+1;
            for (size_t& idex : to_retain) {
                idex -= start_row;
            }
        }
    }

    // Setting up the output structures.
    Rcpp::NumericVector libsizes(cslen);
    const size_t final_nrow=to_retain.size();
    Rcpp::NumericVector outgoing(final_nrow), averaged(final_nrow);

    beachmat::output_param oparam(inmat, false, true);
    oparam.set_chunk_dim(final_nrow, 1); // pure-column chunks for random access, if HDF5.
    auto omat=beachmat::create_numeric_output(final_nrow, cslen, oparam);

    auto lbIt=libsizes.begin();
    size_t out_c=0;
    for (const auto& in_c : csubout) {

        // Extracting the column, subsetting the rows.
        auto inIt=in->get_const_col(in_c, incoming.begin(), start_row, end_row);
        {
            auto oIt=outgoing.begin();
            for (const auto& tr : to_retain) { 
                (*oIt)=*(inIt + tr);
                ++oIt;
            }
        }
           
        // Dividing by the library size. 
        const double& curlib=((*lbIt)=std::accumulate(outgoing.begin(), outgoing.end(), 0.0));
        if (curlib < 0.00000001) {
            throw std::runtime_error("cells should have non-zero library sizes");
        }
        for (double& out : outgoing) {
            out/=curlib;
        }
        ++lbIt;

        omat->set_col(out_c, outgoing.begin());
        ++out_c;

        // Adding to the average.
        {
            auto oIt=outgoing.begin();
            for (auto& ave : averaged) { 
                ave+=(*oIt);
                ++oIt;
            }
        }
    }

    /* Expanding the average vector back to the dimensions spanned by subset_row 
     * (i.e., before removing all-zeroes). This ensures pseudo-cells are comparable 
     * between clusters. Done in C++ to avoid needing to pass back the subset vector.
     */
    Rcpp::IntegerVector kept(to_retain.begin(), to_retain.end());
    return Rcpp::List::create(libsizes, omat->yield(), averaged, kept);
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(matrix);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(matrix);
        return subset_and_divide_internal<Rcpp::IntegerVector>(input.get(), matrix, row_subset, col_subset);
    } else {
        auto input=beachmat::create_numeric_matrix(matrix);
        return subset_and_divide_internal<Rcpp::NumericVector>(input.get(), matrix, row_subset, col_subset);
    }
    END_RCPP
}

/******************** UTILITIES ************************/

int calculate_ratios(const Rcpp::NumericVector& top, const Rcpp::NumericVector& bottom, Rcpp::NumericVector& ratios) {
    int invalid=0;
    auto tIt=top.begin();
    auto bIt=bottom.begin();
    for (auto& r : ratios) {
        if ((*bIt)==0) {
            r=R_PosInf; 
            if ((*tIt)==0) { 
                ++invalid;
            }
        } else {
            r=(*tIt)/(*bIt);
        }
        ++bIt;
        ++tIt;
    }
    return invalid;
}

double compute_median (Rcpp::NumericVector& values, const size_t invalid=0) {
    const size_t validgenes=values.size()-invalid;
    const bool is_even=bool(validgenes%2==0);
    const size_t halfway=size_t(validgenes/2);
    
    // Computing the median (faster than partial sort).
    std::nth_element(values.begin(), values.begin()+halfway, values.end());
    if (is_even) {
        double medtmp=values[halfway];
        std::nth_element(values.begin(), values.begin()+halfway-1, values.end());
        return (medtmp+values[halfway-1])/2;
    } else {
        return values[halfway];
    }       

/*
    std::partial_sort(values, values+halfway+1, values+ngenes)
    if (is_even) {
        return (values[halfway]+values[halfway-1])/2;
    } else 
        return values[halfway];
    } 
*/
}

/*** A function to estimate the pooled size factors and construct the linear equations. ***/


SEXP forge_system (SEXP exprs, SEXP ref, SEXP ordering, SEXP poolsizes) {
    BEGIN_RCPP
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t ngenes=emat->get_nrow();
    const size_t ncells=emat->get_ncol();
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
   
    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(poolsizes);
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        throw std::runtime_error("sizes should be a non-empty integer vector"); 
    }
    int last_size=-1, total_SIZE=0;
    for (const auto& SIZE : pool_sizes) { 
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (SIZE < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
        last_size=SIZE;
    }

    // Checking pseudo cell.
    Rcpp::NumericVector pseudo_cell(ref);
    if (ngenes!=pseudo_cell.size()) { throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); }

    // Checking ordering.
    Rcpp::IntegerVector order(ordering);
    if (order.size() < ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    for (const auto& o : order) { 
        if (o < 0 || o > ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    // Filling up the cell vector.
    Rcpp::NumericVector all_collected(last_size*ngenes);
    std::deque<Rcpp::NumericVector::const_iterator> collected;
    auto acIt=all_collected.begin();
    collected.push_back(acIt); // unfilled first vector, which gets dropped and refilled in the first iteration anyway.
    acIt+=ngenes;

    auto orIt_tail=order.begin();
    Rcpp::NumericVector tmp(emat->get_nrow());
    for (int s=1; s<last_size; ++s, ++orIt_tail, acIt+=ngenes) {
        auto colIt=emat->get_const_col(*orIt_tail, acIt);
        collected.push_back(colIt); 
    }

    // Setting up the output vectors.
    Rcpp::IntegerVector row_num(total_SIZE*ncells), col_num(total_SIZE*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    // Various other bits and pieces.
    Rcpp::NumericVector combined(ngenes), ratios(ngenes);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();

    // Running through the sliding windows.
    for (size_t win=0; win<ncells; ++win, ++orIt) {
        std::fill(combined.begin(), combined.end(), 0);
        
        // Dropping the column that we've moved past, adding the next column.
        if (acIt==all_collected.end()) { 
            acIt=all_collected.begin();
        }
        collected.pop_front();
        auto curIt=emat->get_const_col(*orIt_tail, acIt);
        collected.push_back(curIt);
        ++orIt_tail;
        acIt+=ngenes;

        int index=0;
        int rownum=win; // Setting the row so that all pools with the same SIZE form consecutive equations.
        for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt, rownum+=ncells) { 
            const int& SIZE=(*psIt);
            std::fill(rowIt, rowIt+SIZE, rownum);
            rowIt+=SIZE;
            std::copy(orIt, orIt+SIZE, colIt);
            colIt+=SIZE;

            for (; index<SIZE; ++index) {
                auto ceIt=collected[index];
                for (auto cIt=combined.begin(); cIt!=combined.end(); ++cIt, ++ceIt) {
                    (*cIt)+=(*ceIt);
                }
            }
           
            // Computing the median ratio against the reference.
            int invalid=calculate_ratios(combined, pseudo_cell, ratios);
            pool_factor[rownum]=compute_median(ratios, invalid);
        }

    }    

    return Rcpp::List::create(row_num, col_num, pool_factor);
    END_RCPP
}

/*** Creating the linear system based on the nearest neighbours. ***/

template<class M>
void load_and_average_cache (M mat, Rcpp::NumericVector& cache, std::vector<Rcpp::NumericVector::const_iterator>& cache_iters, 
        Rcpp::NumericVector& pseudo_cell, const size_t i, Rcpp::IntegerMatrix::iterator order_start, const size_t nn) { 

    const size_t ngenes=mat->get_nrow();
    auto cIt=cache.begin();
    cache_iters[0]=mat->get_const_col(i, cIt);
    cIt+=ngenes;

    for (size_t j=0; j<nn; ++j) {
        cache_iters[j+1]=mat->get_const_col(*order_start, cIt);
        ++order_start;
        cIt+=ngenes;
    }   

    // Setting up the pseudo cell.
    std::copy(cache_iters[0], cache_iters[0]+ngenes, pseudo_cell.begin());
    for (size_t j=0; j<nn; ++j) { 
        auto ciIt=cache_iters[j+1];
        auto psIt=pseudo_cell.begin();
        for (size_t i=0; i<ngenes; ++i) { 
            (*psIt)+=(*ciIt);
            ++ciIt;
            ++psIt;
        }
    }

    const double ncells=nn+1;
    for (auto& ps : pseudo_cell) {
        ps/=ncells;
    }
    return;    
} 

SEXP forge_NN_system (SEXP exprs, SEXP ref, SEXP nearest, SEXP poolsizes) {
    BEGIN_RCPP
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t ngenes=emat->get_nrow();
    const size_t ncells=emat->get_ncol();
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
   
    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(poolsizes);
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        throw std::runtime_error("sizes should be a non-empty integer vector"); 
    }
    int last_size=-1, total_SIZE=0;
    for (const auto& SIZE : pool_sizes) { 
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (SIZE < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
        last_size=SIZE;
    }

    // Checking ordering.
    Rcpp::IntegerMatrix order(nearest);
    if (order.ncol() != ncells) { throw std::runtime_error("number of columns in NN matrix should be equal to number of cells"); }
    const size_t nn=order.nrow();
    if (last_size - 1L > nn) { throw std::runtime_error("NN matrix has fewer rows than largest pool"); }
    for (const auto& o : order) { 
        if (o < 0 || o > ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    // Constructing the cache parameters, and also checking the overall pseudo cell.
    Rcpp::IntegerVector chosen(ref);
    if (chosen.size()!=1) { throw std::runtime_error("pseudo-cell specification should be an integer scalar"); }
    Rcpp::NumericVector exprs_cache(ngenes*(nn+1));
    std::vector<Rcpp::NumericVector::const_iterator> cache_iters(nn+1);
    Rcpp::NumericVector pseudo_cell(ngenes);
    load_and_average_cache(emat.get(), exprs_cache, cache_iters, pseudo_cell, 
            chosen[0], order.begin() + chosen[0]*nn, nn);
   
    // Setting up the output vectors.
    Rcpp::IntegerVector row_num(total_SIZE*ncells), col_num(total_SIZE*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);
    Rcpp::NumericVector combined(ngenes), ratios(ngenes), current_pseudo_cell(ngenes);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();

    // Running through each cell and constructing the pools based on the nearest neighbours.
    // Each pool is normalized to the current pseudo cell, and then the results are scaled
    // based on the normalization of the current pseudo cell to the overall pseudo cell.
    for (size_t cell=0; cell<ncells; ++cell) {
        load_and_average_cache(emat.get(), exprs_cache, cache_iters, 
                current_pseudo_cell, cell, orIt, nn);

        // Computing the ratio of the two current and overall pseudocell.
        int invalid=calculate_ratios(current_pseudo_cell, pseudo_cell, ratios);
        const double scaling=compute_median(ratios, invalid);

        // Copying self first before iteratively adding the nearest neighbours.
        std::copy(cache_iters[0], cache_iters[0]+ngenes, combined.begin());
        int index=1;
        int rownum=cell; 

        for (const auto& SIZE : pool_sizes) { 
            std::fill(rowIt, rowIt+SIZE, rownum);
            rowIt+=SIZE;
            (*colIt)=cell;
            std::copy(orIt, orIt+SIZE-1, colIt+1);
            colIt+=SIZE;

            for (; index<SIZE; ++index) {
                auto ceIt=cache_iters[index];
                for (auto& co : combined) {
                    co+=(*ceIt);
                    ++ceIt;
                }
            }
           
            // Computing the ratio against the current reference.
            int invalid=calculate_ratios(combined, current_pseudo_cell, ratios);
            pool_factor[rownum]=compute_median(ratios, invalid) * scaling;

            // Bumping the row number so that all pools with the same SIZE form consecutive equations.
            rownum+=ncells;
        }
        orIt+=nn;
    }    

    return Rcpp::List::create(row_num, col_num, pool_factor, pseudo_cell);
    END_RCPP
}


