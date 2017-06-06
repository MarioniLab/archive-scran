#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

#include "scran.h"

/*** A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 *** The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 *** computed as colSums(MAT[row_subset,col_subset]).
 ***/

template <class V, class INMAT>
SEXP subset_and_divide_internal(const INMAT in, SEXP inmat, SEXP row_subset, SEXP col_subset) {
    // Checking row subset vector
    subset_values rsubout=check_subset_vector(row_subset, in->get_nrow());
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    // Checking column subset vector
    subset_values csubout=check_subset_vector(col_subset, in->get_ncol());
    const int cslen=csubout.first;
    const int* csptr=csubout.second;

    // Setting up the output structures.
    auto out=beachmat::create_numeric_output(rslen, cslen, inmat, false, true);
    Rcpp::NumericVector libsizes(cslen);
    V incoming(in->get_nrow());
    Rcpp::NumericVector outgoing(rslen);

    auto lIt=libsizes.begin();
    for (int cs=0; cs<cslen; ++cs, ++lIt) {

        in->get_col(csptr[cs], incoming.begin());
        auto oIt=outgoing.begin();
        for (int rs=0; rs<rslen; ++rs, ++oIt) {
            (*oIt)=*(incoming.begin() + rsptr[rs]);
        }
            
        const double& curlib=(*lIt=std::accumulate(outgoing.begin(), outgoing.end(), 0.0));
        if (curlib < 0.00000001) {
            throw std::runtime_error("cells should have non-zero library sizes");
        }

        for (auto oIt=outgoing.begin(); oIt!=outgoing.end(); ++oIt) { 
            (*oIt)/=curlib;
        }

        out->fill_col(cs, outgoing.begin());
    }

    Rcpp::List output(2);
    output[0] = libsizes;
    output[1] = out->yield();
    return SEXP(output);
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

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

SEXP forge_system (SEXP exprs, SEXP ordering, SEXP sizes, SEXP ref) {
    BEGIN_RCPP

    // Checking input matrix.
    auto emat=beachmat::create_numeric_matrix(exprs);
    const size_t ncells=emat->get_ncol();
    const size_t ngenes=emat->get_nrow();
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }

    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(sizes);
    const size_t nsizes=pool_sizes.size();
    if (!nsizes) {
        throw std::runtime_error("sizes should be a non-empty integer vector"); 
    }
    int last_size=-1;
    size_t total_SIZE=0;
    for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt) { 
        const int& SIZE=*psIt;
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (SIZE < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
        last_size=SIZE;
    }

    // Checking reference and ordering.
    Rcpp::NumericVector pseudo_cell(ref);
    if (ngenes!=pseudo_cell.size()) { throw std::runtime_error("length of reference vector is inconsistent with number of cells"); }
    Rcpp::IntegerVector order(ordering);
    if (order.size() < ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }

    // Setting up the output matrix.
    Rcpp::IntegerVector row_num(total_SIZE*ncells), col_num(total_SIZE*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    // Various bits and pieces.
    std::vector<double> combined(ngenes), ratios(ngenes);
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();
   
    // Filling up the cell vector.
    Rcpp::NumericVector all_collected(last_size*ngenes);
    std::deque<Rcpp::NumericVector::iterator> collected;
    auto acIt=all_collected.begin();
    collected.push_back(acIt);
    acIt+=ngenes;

    auto orIt_tail=order.begin();
    for (int s=0; s<last_size-1; ++s, ++orIt_tail, acIt+=ngenes) {
        emat->get_col(*orIt_tail, acIt);
        collected.push_back(acIt);
    }

    // Running through the sliding windows.
    for (size_t win=0; win<ncells; ++win, ++orIt) {
        std::fill(combined.begin(), combined.end(), 0);
        
        // Dropping the first element, adding the last element.
        auto tmp=collected.front();
        collected.pop_front();
        emat->get_col(*orIt_tail, tmp);
        collected.push_back(tmp);
        ++orIt_tail;

        int index=0;
        int rownum=win; // Setting the row so that all rows with the same SIZE are consecutive.
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
           
            // Computing the ratio against the reference.
            auto rIt=ratios.begin(), cIt=combined.begin();
            for (auto pcIt=pseudo_cell.begin(); pcIt!=pseudo_cell.end(); ++pcIt, ++rIt, ++cIt) {
                (*rIt)=(*cIt)/(*pcIt);
            }

            // Computing the median (faster than partial sort).
            std::nth_element(ratios.begin(), ratios.begin()+halfway, ratios.end());
            if (is_even) {
                double medtmp=ratios[halfway];
                std::nth_element(ratios.begin(), ratios.begin()+halfway-1, ratios.end());
                pool_factor[rownum]=(medtmp+ratios[halfway-1])/2;
            } else {
                pool_factor[rownum]=ratios[halfway];
            }       
        }

/*
        std::partial_sort(combined, combined+halfway+1, combined+ngenes);
        if (is_even) {
            ofptr[cell]=(combined[halfway]+combined[halfway-1])/2;
        } else {
            ofptr[cell]=combined[halfway];
        } 
*/
    }    

    Rcpp::List output(3);
    output[0]=row_num;
    output[1]=col_num;
    output[2]=pool_factor;
    return SEXP(output);

    END_RCPP
}

