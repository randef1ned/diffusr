// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/diffusr.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// heat_diffusion_
Eigen::MatrixXd heat_diffusion_(const Eigen::MatrixXd& v0, const Eigen::MatrixXd& W, const double t);
static SEXP _diffusr_heat_diffusion__try(SEXP v0SEXP, SEXP WSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(heat_diffusion_(v0, W, t));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _diffusr_heat_diffusion_(SEXP v0SEXP, SEXP WSEXP, SEXP tSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_diffusr_heat_diffusion__try(v0SEXP, WSEXP, tSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// stoch_col_norm_
Eigen::MatrixXd stoch_col_norm_(const Eigen::MatrixXd& W);
static SEXP _diffusr_stoch_col_norm__try(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(stoch_col_norm_(W));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _diffusr_stoch_col_norm_(SEXP WSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_diffusr_stoch_col_norm__try(WSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// laplacian_
Eigen::MatrixXd laplacian_(const Eigen::MatrixXd& W);
static SEXP _diffusr_laplacian__try(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(laplacian_(W));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _diffusr_laplacian_(SEXP WSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_diffusr_laplacian__try(WSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mrwr_
Eigen::MatrixXd mrwr_(const Eigen::MatrixXd& p0, const Eigen::MatrixXd& W, const double r, const double thresh, const int niter, const bool do_analytical);
static SEXP _diffusr_mrwr__try(SEXP p0SEXP, SEXP WSEXP, SEXP rSEXP, SEXP threshSEXP, SEXP niterSEXP, SEXP do_analyticalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_analytical(do_analyticalSEXP);
    rcpp_result_gen = Rcpp::wrap(mrwr_(p0, W, r, thresh, niter, do_analytical));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _diffusr_mrwr_(SEXP p0SEXP, SEXP WSEXP, SEXP rSEXP, SEXP threshSEXP, SEXP niterSEXP, SEXP do_analyticalSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_diffusr_mrwr__try(p0SEXP, WSEXP, rSEXP, threshSEXP, niterSEXP, do_analyticalSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// neighbors_
Rcpp::List neighbors_(const Rcpp::IntegerVector& node_idxs, const Rcpp::NumericMatrix& W, const int k);
static SEXP _diffusr_neighbors__try(SEXP node_idxsSEXP, SEXP WSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type node_idxs(node_idxsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(neighbors_(node_idxs, W, k));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _diffusr_neighbors_(SEXP node_idxsSEXP, SEXP WSEXP, SEXP kSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_diffusr_neighbors__try(node_idxsSEXP, WSEXP, kSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _diffusr_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Eigen::MatrixXd(*heat_diffusion_)(const Eigen::MatrixXd&,const Eigen::MatrixXd&,const double)");
        signatures.insert("Eigen::MatrixXd(*stoch_col_norm_)(const Eigen::MatrixXd&)");
        signatures.insert("Eigen::MatrixXd(*laplacian_)(const Eigen::MatrixXd&)");
        signatures.insert("Eigen::MatrixXd(*mrwr_)(const Eigen::MatrixXd&,const Eigen::MatrixXd&,const double,const double,const int,const bool)");
        signatures.insert("Rcpp::List(*neighbors_)(const Rcpp::IntegerVector&,const Rcpp::NumericMatrix&,const int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _diffusr_RcppExport_registerCCallable() { 
    R_RegisterCCallable("diffusr", "_diffusr_heat_diffusion_", (DL_FUNC)_diffusr_heat_diffusion__try);
    R_RegisterCCallable("diffusr", "_diffusr_stoch_col_norm_", (DL_FUNC)_diffusr_stoch_col_norm__try);
    R_RegisterCCallable("diffusr", "_diffusr_laplacian_", (DL_FUNC)_diffusr_laplacian__try);
    R_RegisterCCallable("diffusr", "_diffusr_mrwr_", (DL_FUNC)_diffusr_mrwr__try);
    R_RegisterCCallable("diffusr", "_diffusr_neighbors_", (DL_FUNC)_diffusr_neighbors__try);
    R_RegisterCCallable("diffusr", "_diffusr_RcppExport_validate", (DL_FUNC)_diffusr_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_diffusr_heat_diffusion_", (DL_FUNC) &_diffusr_heat_diffusion_, 3},
    {"_diffusr_stoch_col_norm_", (DL_FUNC) &_diffusr_stoch_col_norm_, 1},
    {"_diffusr_laplacian_", (DL_FUNC) &_diffusr_laplacian_, 1},
    {"_diffusr_mrwr_", (DL_FUNC) &_diffusr_mrwr_, 6},
    {"_diffusr_neighbors_", (DL_FUNC) &_diffusr_neighbors_, 3},
    {"_diffusr_RcppExport_registerCCallable", (DL_FUNC) &_diffusr_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_diffusr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
