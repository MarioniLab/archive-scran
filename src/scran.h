#ifndef SCRAN_H
#define SCRAN_H

#include "R.h"
#include "Rinternals.h"
#include <stdexcept>
#include <algorithm>

extern "C" {

SEXP forge_system (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP shuffle_scores (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#include "utils.hpp"

#endif 
