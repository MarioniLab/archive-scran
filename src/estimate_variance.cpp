#include "scran.h"
#include "run_dormqr.h"

SEXP estimate_variance (SEXP qr, SEXP qraux, SEXP exprs, SEXP subset) {
    BEGIN_RCPP

    // Setting up for Q-based multiplication.
    run_dormqr multQ(qr, qraux, 'T');
    const int ncoefs=multQ.get_ncoefs();
    const int ncells=multQ.get_nobs();

    // Setting up expression data.
    auto emat=beachmat::create_numeric_matrix(exprs);
    if (ncells!=int(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    }
    auto subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.size();

    // Setting up output objects.
    Rcpp::NumericVector means(slen), vars(slen);
    auto mIt=means.begin(), vIt=vars.begin();
    Rcpp::NumericVector tmp(ncells);

    // Running through each gene and reporting its variance and mean.
    for (auto sIt=subout.begin(); sIt!=subout.end(); ++sIt, ++mIt, ++vIt) {
        emat->get_row(*sIt, tmp.begin());
        (*mIt)=std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells;

        std::copy(tmp.begin(), tmp.end(), multQ.rhs.begin());
        multQ.run();

        double& curvar=(*vIt);
        for (int c=ncoefs; c<ncells; ++c) { // only using the residual effects.
            curvar += multQ.rhs[c]*multQ.rhs[c];
        }
        curvar /= ncells - ncoefs;
    }

    return Rcpp::List::create(means, vars);
    END_RCPP
}
