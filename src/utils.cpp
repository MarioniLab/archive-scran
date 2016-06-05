#include "scran.h"

/* Various utility functions that don't really warrant their own page, 
 * but are helpful for avoiding construction of large temporaries. */

typedef std::pair<const int, const int*> subset_values;
subset_values check_subset_vector(SEXP subset, int maxdim) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }
    const int slen=LENGTH(subset);
    const int* sptr=INTEGER(subset);
    for (int s=0; s<slen; ++s) {
        if (sptr[s] < 0 || sptr[s] >= maxdim) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return std::make_pair(slen, sptr);
}

/* A function to (a) subset by row, (b) subset by column, and (c) divide through by the library sizes. 
 * The output is equivalent to t(t(MAT[row_subset,col_subset])/lib.sizes) where lib.sizes itself is
 * computed as colSums(MAT[row_subset,col_subset]).
 */

template <typename T>
SEXP subset_and_divide_internal(const T* ptr, const matrix_info& MAT, SEXP row_subset, SEXP col_subset) {
    // Checking row subset vector
    subset_values rsubout=check_subset_vector(row_subset, int(MAT.nrow));
    const int rslen=rsubout.first;
    const int* rsptr=rsubout.second;

    // Checking column subset vector
    subset_values csubout=check_subset_vector(col_subset, int(MAT.ncol));
    const int cslen=csubout.first;
    const int* csptr=csubout.second;

    SEXP output=PROTECT(allocVector(VECSXP, 2));
    try {
        SET_VECTOR_ELT(output, 0, allocVector(REALSXP, cslen));
        double* olptr=REAL(VECTOR_ELT(output, 0));

        SET_VECTOR_ELT(output, 1, allocMatrix(REALSXP, rslen, cslen));
        double* onptr=REAL(VECTOR_ELT(output, 1));
        const T* curptr;

        for (int cs=0; cs<cslen; ++cs) {
            curptr=ptr + MAT.nrow*csptr[cs];

            double& curlib=(olptr[cs]=0);
            for (int rs=0; rs<rslen; ++rs) { 
                curlib+=curptr[rsptr[rs]];
            }
            if (curlib < 0.00000001) {
                throw std::runtime_error("cells should have non-zero library sizes");
            }

            for (int rs=0; rs<rslen; ++rs) { 
                onptr[rs] = curptr[rsptr[rs]]/curlib;
            }
            onptr+=rslen;
        }

    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP subset_and_divide(SEXP matrix, SEXP row_subset, SEXP col_subset) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return subset_and_divide_internal<int>(MAT.iptr, MAT, row_subset, col_subset);
    } else {
        return subset_and_divide_internal<double>(MAT.dptr, MAT, row_subset, col_subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}


