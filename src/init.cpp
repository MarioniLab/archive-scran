#include "scran.h"
#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(forge_system, 6),
    REGISTER(shuffle_scores, 8),
    REGISTER(get_null_rho, 2),
    REGISTER(get_null_rho_design, 4),
    REGISTER(compute_rho, 4),
    REGISTER(auto_shuffle, 2),
    REGISTER(compute_cordist, 2),
    REGISTER(subset_and_divide, 3),
    {NULL, NULL, 0}
};

void attribute_visible R_init_scran(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

