#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "Rcpp.h"

#include "scran.h"

SEXP estimate_variance (SEXP qr, SEXP qraux, SEXP exprs, SEXP subset) {
    BEGIN_RCPP

    // Setting up for Q-based multiplication.
    Rcpp::NumericMatrix QR(qr);
    Rcpp::NumericVector QRaux(qraux);
    if (QRaux.size()!=QR.ncol()) {
        throw std::runtime_error("QR auxiliary vector should be of length 'ncol(Q)'");
    }
    const double * qrptr=NULL, * qrxptr=NULL;
    if (QR.size()) { 
        qrptr=&(QR[0]);
    }
    if (QRaux.size()) { 
       qrxptr=&(QRaux[0]);
    }
    const int ncoefs=QR.ncol();
    const int ncells=QR.nrow();
    run_dormqr multQ(ncells, ncoefs, qrptr, qrxptr, 'T');

    // Setting up expression data.
    auto emat=beachmat::create_numeric_matrix(exprs);
    if (ncells!=int(emat->get_ncol())) {
        throw std::runtime_error("number of rows of QR matrix not equal to number of cells");
    }
    subset_values subout=check_subset_vector(subset, emat->get_nrow());
    const size_t slen=subout.first;
    const int* sptr=subout.second;

    // Setting up output objects.
    Rcpp::NumericVector means(slen), vars(slen);
    auto mIt=means.begin(), vIt=vars.begin();
    Rcpp::NumericVector tmp(ncells);

    // Running through each gene and reporting its variance and mean.
    for (int s=0; s<slen; ++s, ++mIt, ++vIt) {
        emat->get_row(sptr[s], tmp.begin());
        (*mIt)=std::accumulate(tmp.begin(), tmp.end(), 0.0)/ncells;

        std::copy(tmp.begin(), tmp.end(), multQ.rhs);
        multQ.run();

        double& curvar=(*vIt);
        for (int c=ncoefs; c<ncells; ++c) { // only using the residual effects.
            curvar += multQ.rhs[c]*multQ.rhs[c];
        }
        curvar /= ncells - ncoefs;
    }

    Rcpp::List output(2);
    output[0]=means;
    output[1]=vars;
    return output;
    END_RCPP
}
