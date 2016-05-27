#include "scran.h"

double rho_mult (double Ncells) {
    return 6/(Ncells*(Ncells*Ncells-1));
}

/*** Null distribution estimation without a design matrix. ***/
SEXP get_null_rho (SEXP cells, SEXP iters) try {
    if (!isInteger(cells) || LENGTH(cells)!=1)  {
        throw std::runtime_error("number of cells should be an integer scalar"); 
    }
    const int Ncells=asInteger(cells);
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }
    if (!isInteger(iters) || LENGTH(iters)!=1)  {
        throw std::runtime_error("number of iterations should be an integer scalar"); 
    }
    const int Niters=asInteger(iters);
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }

    int* rankings=(int*)R_alloc(Ncells, sizeof(int));
    int cell;
    for (cell=0; cell<Ncells; ++cell) { rankings[cell]=cell; }    

    SEXP output=PROTECT(allocVector(REALSXP, Niters));
    try {
        double* optr=REAL(output);
        Rx_random_seed myseed;
        double tmp, tmpdiff;
        const double mult=rho_mult(Ncells);

        for (int it=0; it<Niters; ++it) {
            Rx_shuffle(rankings, rankings + Ncells);
            tmp=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmpdiff=rankings[cell]-cell;
                tmp+=tmpdiff*tmpdiff;
            }
            tmp*=mult;
            optr[it]=1-tmp;            
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

/*** Null distribution estimation with a design matrix. ***/
SEXP get_null_rho_design(SEXP design, SEXP coef, SEXP obs, SEXP iters) try {
    if (!isInteger(coef) || LENGTH(coef)!=1)  {
        throw std::runtime_error("number of coef should be an integer scalar"); 
    }
    const int Ncoef=asInteger(coef);
    if (!isInteger(obs) || LENGTH(obs)!=1)  {
        throw std::runtime_error("number of obs should be an integer scalar"); 
    }
    const int Nobs=asInteger(obs);
    if (!isInteger(iters) || LENGTH(iters)!=1)  {
        throw std::runtime_error("number of iterations should be an integer scalar"); 
    }
    const int Niters=asInteger(iters);
    if (Niters <= 0) { throw std::runtime_error("number of iterations should be positive"); }

    // Setting up the Q matrix (in compressed form + auxiliaries).
    if (!isReal(design) || LENGTH(design)!=Nobs*Ncoef) {
        throw std::runtime_error("design matrix dimensions are not consistent with number of observations");
    }
    const double* xptr=REAL(design);
    double* qptr=(double*)R_alloc(LENGTH(design), sizeof(double));
    double* qxptr=(double*)R_alloc(Ncoef, sizeof(double)); /* number of elementary reflectors is the minimum of Nobs, Ncoef */
    {
        // Getting the QR decomposition after a workspace query.
        for (int i=0; i<LENGTH(design); ++i) { qptr[i]=xptr[i]; }
        int info, lwork=-1;
        double tmpwork=0;
        F77_CALL(dgeqrf)(&Nobs, &Ncoef, qptr, &Nobs, qxptr, &tmpwork, &lwork, &info);

        lwork=int(tmpwork + 0.5);
        double* work=(double*)R_alloc(lwork, sizeof(double));
        F77_CALL(dgeqrf)(&Nobs, &Ncoef, qptr, &Nobs, qxptr, work, &lwork, &info);
    }

    // Workspace query for 'dormqr'.
    double* effects=(double*)R_alloc(Nobs, sizeof(double));
    const char side='L', trans='N';
    const int nCols=1;
    int info, lwork=-1;
    double* work;
    {
        double tmpwork=0;
        F77_CALL(dormqr)(&side, &trans, &Nobs, &nCols, &Ncoef,
                qptr, &Nobs, qxptr, effects, &Nobs,
                &tmpwork, &lwork, &info); 
        if (info) { 
            throw std::runtime_error("workspace query failed for 'dormqr'");
        }
        lwork=int(tmpwork+0.5);
        work=(double*)R_alloc(lwork, sizeof(double));
    }

    SEXP output=PROTECT(allocVector(REALSXP, Niters));
    try {
        double* optr=REAL(output);
        Rx_random_seed myseed;
    
        std::deque<std::pair<double, int> > collected1(Nobs), collected2(Nobs);
        std::deque<int> rank1(Nobs), rank2(Nobs);
        const double mult=rho_mult(Nobs);

        // Simulating residuals, using the Q-matrix to do it.
        // We set the main effects to zero (hence, starting from "Ncoefs") and simulate normals for the residual effects.
        // We then use this to reconstruct the residuals themselves, and then compute correlations between them.
        int mode, row, col;
        double tmpdiff;

        for (int it=0; it<Niters; ++it) {
            for (mode=0; mode<2; ++mode) {
                for (col=0; col<Ncoef; ++col) {
                    effects[col]=0;
                }
                for (col=Ncoef; col<Nobs; ++col) {
                    effects[col]=norm_rand();
                }

                // Computing the residuals.
                F77_CALL(dormqr)(&side, &trans, &Nobs, &nCols, &Ncoef,
                        qptr, &Nobs, qxptr, effects, &Nobs,
                        work, &lwork, &info); 
                if (info) { 
                    throw std::runtime_error("residual calculations failed for 'dormqr'");
                }

                // Sorting.
                std::deque<std::pair<double, int> >& current=(mode ? collected1 : collected2);
                for (row=0; row<Nobs; ++row) {
                    current[row].first=effects[row];
                    current[row].second=row;
                }
                std::sort(current.begin(), current.end());
                std::deque<int>& rank=(mode ? rank1 : rank2);
                for (row=0; row<Nobs; ++row) {
                    rank[current[row].second]=row;
                }
            }

            // Computing the squared difference in the ranks.
            double& rho=(optr[it]=0);
            for (row=0; row<Nobs; ++row) {
                tmpdiff=rank1[row]-rank2[row];
                rho+=tmpdiff*tmpdiff;
            }
            rho*=mult;
            rho=1-rho;
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

/*** Estimating correlations (without expanding into a matrix to do so via 'cor'). ***/
SEXP compute_rho(SEXP g1, SEXP g2, SEXP cells, SEXP rankings) try {
    if (!isInteger(cells) || LENGTH(cells)!=1)  {
        throw std::runtime_error("number of cells should be an integer scalar"); 
    }
    const int Ncells=asInteger(cells);
    if (Ncells <= 1) { throw std::runtime_error("number of cells should be greater than 2"); }

    if (!isInteger(g1) || !isInteger(g2)) { 
        throw std::runtime_error("gene indices must be integer vectors");
    }
    const int Npairs=LENGTH(g1);
    if (Npairs!=LENGTH(g2)) { 
        throw std::runtime_error("gene index vectors must be of the same length"); 
    }
    const int *g1ptr=INTEGER(g1), *g2ptr=INTEGER(g2);
    
    if (!isInteger(rankings)) { 
        throw std::runtime_error("ranking matrix must be integer");
    }
    const int Ngenes=LENGTH(rankings)/Ncells;
    if (Ngenes*Ncells!=LENGTH(rankings)) { 
        throw std::runtime_error("dimensions of ranking matrix are not consistent with number of cells"); 
    }
    const int* rptr=INTEGER(rankings);

    SEXP output=PROTECT(allocVector(REALSXP, Npairs));
    try {
        double* orptr=REAL(output);
        
        int off1, off2, cell;
        double tmp, tmp2;
        const double mult=rho_mult(Ncells); 

        for (int p=0; p<Npairs; ++p) {
            off1=(g1ptr[p]-1)*Ncells;
            off2=(g2ptr[p]-1)*Ncells;
            if (g1ptr[p] < 1 || g1ptr[p] > Ngenes) {
                throw std::runtime_error("first gene index is out of range");
            }
            if (g2ptr[p] < 1 || g2ptr[p] > Ngenes) {
                throw std::runtime_error("second gene index is out of range");
            }

            // Computing the correlation.
            tmp=0;
            for (cell=0; cell<Ncells; ++cell) {
                tmp2=rptr[off1+cell] - rptr[off2+cell];
                tmp+=tmp2*tmp2;
            }
            tmp*=mult;
            orptr[p]=1-tmp;
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


