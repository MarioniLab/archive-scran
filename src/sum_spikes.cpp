#include "scran.h"

template<class M>
SEXP sum_spikes_internal(M mat, Rcpp::IntegerVector spikedex) {
    const size& ngenes = mat->nrow();
    const size& ncells = mat->ncol();
    int last=-1;
    for (const auto& curdex : spikedex) { 
        if (curdex <= last || curdex<0 || curdex>=ngenes) {
            throw std::runtime_error("'spikedex' should contain sorted indices in [0, ngenes)");
        } 
    }

    Rcpp::NumericVector incoming(ngenes);
    Rcpp::NumericVector output(ncells);
    auto oIt=output.begin();
    for (size_t c=0; c<ncells; ++c, ++oIt) {
        auto iIt=mat->get_const_col(c, incoming.begin());
        for (const auto& curdex : spikedex) {
            (*oIt)+=*(iIt + curdex);
        }
    }

    return output;
}

SEXP sum_spikes(SEXP mat, SEXP spikedex) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(exprs);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return sum_spikes_internal(mat, spikedex);
    } else {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return sum_spikes_internal(mat, spikedex);
    }
    END_RCPP
}
