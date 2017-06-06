#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

#include "scran.h"

/*** A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 *** The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 *** computed as colSums(MAT[row_subset,col_subset]).
 ***/

template <class V, class M>
SEXP subset_and_divide_internal(const M in, SEXP row_subset, SEXP col_subset) {
    // Checking row subset vector
    subset_values rsubout=check_subset_vector(row_subset, in->get_nrow());
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    // Checking column subset vector
    subset_values csubout=check_subset_vector(col_subset, in->get_ncol());
    const int cslen=csubout.first;
    const int* csptr=csubout.second;

    // Setting up the output structures.
    Rcpp::NumericVector libsizes(cslen);
    V incoming(in->get_nrow());
    Rcpp::NumericVector outgoing(rslen);
    Rcpp::NumericVector pseudocell(rslen);

    for (int cs=0; cs<cslen; ++cs) {

        in->get_col(csptr[cs], incoming.begin());
        auto oIt=outgoing.begin();
        for (int rs=0; rs<rslen; ++rs, ++oIt) {
            (*oIt)=*(incoming.begin() + rsptr[rs]);
        }
            
        const double& curlib=(libsizes[cs]=std::accumulate(outgoing.begin(), outgoing.end(), 0.0));
        if (curlib < 0.00000001) {
            throw std::runtime_error("cells should have non-zero library sizes");
        }

        auto pcIt=pseudocell.begin();
        for (auto oIt=outgoing.begin(); oIt!=outgoing.end(); ++oIt, ++pcIt) { 
            (*pcIt)+=(*oIt)/curlib;
        }
    }

    // Taking the average to form the pseudocell.
    for (auto pcIt=pseudocell.begin(); pcIt!=pseudocell.end(); ++pcIt) {
        (*pcIt)/=cslen;
    }

    Rcpp::List output(2);
    output[0] = libsizes;
    output[1] = pseudocell;
    return SEXP(output);
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(matrix);
    if (rtype==INTSXP) {
        auto input=beachmat::create_integer_matrix(matrix);
        return subset_and_divide_internal<Rcpp::IntegerVector>(input.get(), row_subset, col_subset);
    } else {
        auto input=beachmat::create_numeric_matrix(matrix);
        return subset_and_divide_internal<Rcpp::NumericVector>(input.get(), row_subset, col_subset);
    }
    END_RCPP
}

/*** A function to estimate the pooled size factors and construct the linear equations. ***/

template<class M>
SEXP forge_system_internal (const M emat, SEXP subset_row, SEXP subset_col, SEXP ref, SEXP libsizes, SEXP ordering, SEXP poolsizes) {
    // Checking row subset vector
    subset_values rsubout=check_subset_vector(subset_row, emat->get_nrow());
    const size_t ngenes=rsubout.first;
    const int* rsptr=rsubout.second;

    subset_values csubout=check_subset_vector(subset_col, emat->get_ncol());
    const size_t ncells=csubout.first;
    const int* csptr=csubout.second;
    if (ncells==0) { throw std::runtime_error("at least one cell required for normalization"); }
   
    // Checking the input sizes.
    Rcpp::IntegerVector pool_sizes(poolsizes);
    const size_t nsizes=pool_sizes.size();
    if (nsizes==0) {
        throw std::runtime_error("sizes should be a non-empty integer vector"); 
    }
    int last_size=-1, total_SIZE=0;
    for (auto psIt=pool_sizes.begin(); psIt!=pool_sizes.end(); ++psIt) { 
        const int& SIZE=*psIt;
        if (SIZE < 1 || SIZE > ncells) { throw std::runtime_error("each element of sizes should be within [1, number of cells]"); }
        if (SIZE < last_size) { throw std::runtime_error("sizes should be sorted"); }
        total_SIZE+=SIZE;
        last_size=SIZE;
    }

    // Checking pseudo cell and library sizes.
    Rcpp::NumericVector pseudo_cell(ref);
    if (ngenes!=pseudo_cell.size()) { throw std::runtime_error("length of pseudo-cell vector is not the same as the number of cells"); }
    Rcpp::NumericVector lib_sizes(libsizes);
    if (ncells!=lib_sizes.size()) { throw std::runtime_error("length of library size vector is not the same as the number of cells"); }

    // Checking ordering.
    Rcpp::IntegerVector order(ordering);
    if (order.size() < ncells*2-1)  { throw std::runtime_error("ordering vector is too short for number of cells"); }
    for (auto oIt=order.begin(); oIt!=order.end(); ++oIt) {
        if (*oIt < 0 || *oIt > ncells) { 
            throw std::runtime_error("elements of ordering vector are out of range");
        }
    }

    // Filling up the cell vector.
    Rcpp::NumericVector all_collected(last_size*ngenes);
    std::deque<Rcpp::NumericVector::iterator> collected;
    auto acIt=all_collected.begin();
    collected.push_back(acIt); // unfilled first vector, which gets dropped and refilled in the first iteration anyway.
    acIt+=ngenes;

    auto orIt_tail=order.begin();
    Rcpp::NumericVector tmp(emat->get_nrow());
    for (int s=0; s<last_size-1; ++s, ++orIt_tail) {
        emat->get_col(csptr[*orIt_tail], tmp.begin());
        collected.push_back(acIt); 
        const double& curlibsize=lib_sizes[*orIt_tail];
        for (size_t g=0; g<ngenes; ++g, ++acIt) {
            (*acIt)=tmp[rsptr[g]]/curlibsize;
        }
    }

    // Setting up the output vectors.
    Rcpp::IntegerVector row_num(total_SIZE*ncells), col_num(total_SIZE*ncells);
    Rcpp::NumericVector pool_factor(nsizes*ncells);

    // Various other bits and pieces.
    std::vector<double> combined(ngenes), ratios(ngenes);
    const bool is_even=bool(ngenes%2==0);
    const int halfway=int(ngenes/2);
    auto rowIt=row_num.begin(), colIt=col_num.begin();
    auto orIt=order.begin();

    // Running through the sliding windows.
    for (size_t win=0; win<ncells; ++win, ++orIt) {
        std::fill(combined.begin(), combined.end(), 0);
        
        // Dropping the column that we've moved past, adding the next column.
        auto acIt=collected.front();
        collected.pop_front();
        collected.push_back(acIt);
        emat->get_col(csptr[*orIt_tail], tmp.begin());
        const double& curlibsize=lib_sizes[*orIt_tail];
        for (size_t g=0; g<ngenes; ++g, ++acIt) {
            (*acIt)=tmp[rsptr[g]]/curlibsize;
        }
        ++orIt_tail;

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
}

SEXP forge_system(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP ref, SEXP libsizes, SEXP ordering, SEXP poolsizes) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto emat=beachmat::create_integer_matrix(exprs);
        return forge_system_internal(emat.get(), subset_row, subset_col, ref, libsizes, ordering, poolsizes);
    } else {
        auto emat=beachmat::create_numeric_matrix(exprs);
        return forge_system_internal(emat.get(), subset_row, subset_col, ref, libsizes, ordering, poolsizes);
    }
    END_RCPP
}

