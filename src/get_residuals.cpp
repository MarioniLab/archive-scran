#include "scran.h"

/* This function computes residuals in a nice and quick manner.
 * It takes a pre-computed QR matrix from qr(), and optionally
 * a subset of integer indices, and returns a matrix of residuals.
 */

SEXP get_residuals(SEXP exprs, SEXP qr, SEXP qraux, SEXP subset) try {
    const matrix_info emat=check_matrix(exprs);
    if (emat.is_integer) {
        throw std::runtime_error("expression matrix must be double-precision"); 
    }
    const double** eptrs=(const double**)R_alloc(emat.ncol, sizeof(const double*));
    if (emat.ncol) {
        eptrs[0]=emat.dptr;
        for (size_t c=1; c<emat.ncol; ++c) {
            eptrs[c]=eptrs[c-1]+emat.nrow;
        }
    }

    // Checking the subset vector.
    subset_values subout=check_subset_vector(subset, emat.nrow);
    const int slen=subout.first;
    const int* sptr=subout.second;
    
    // Checking the QR matrix.
    const matrix_info QR=check_matrix(qr);
    if (QR.is_integer){ 
        throw std::runtime_error("QR matrix must be double-precision");
    }
    if (!isReal(qraux) || size_t(LENGTH(qraux))!=QR.ncol) {
        throw std::runtime_error("QR auxiliary vector should be double-precision and of length 'ncol(Q)'");
    }
    const double* qrxptr=REAL(qraux);
    run_dormqr multQ1(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'T');
    run_dormqr multQ2(QR.nrow, QR.ncol, QR.dptr, qrxptr, 'N');
    
    SEXP output=PROTECT(allocMatrix(REALSXP, slen, emat.ncol));
    try {
        double** optrs=(double**)R_alloc(emat.ncol, sizeof(double*));
        if (emat.ncol) {
            optrs[0]=REAL(output);
            for (size_t c=1; c<emat.ncol; ++c) {
                optrs[c]=optrs[c-1]+slen;
            }
        }

        size_t c;
        double* temporary=(double*)R_alloc(emat.ncol, sizeof(double));
        for (int s=0; s<slen; ++s) {
            for (c=0; c<emat.ncol; ++c) {
                temporary[c]=eptrs[c][sptr[s]];
            }

            multQ1.run(temporary); // Getting main+residual effects.
            for (c=0; c<QR.ncol; ++c) {
                temporary[c]=0; // setting main effects to zero.
            }

            multQ2.run(temporary); // Getting residuals.
            for (c=0; c<emat.ncol; ++c) {
                optrs[c][s]=temporary[c];
            }
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
} catch (std::exception& e) {
    return mkString(e.what());
}

