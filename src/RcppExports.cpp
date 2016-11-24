// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/diffusr.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// heat_diffusion_
Eigen::VectorXd heat_diffusion_(const Eigen::VectorXd& v0, const Eigen::MatrixXd& W, const double b);
static SEXP diffusr_heat_diffusion__try(SEXP v0SEXP, SEXP WSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(heat_diffusion_(v0, W, b));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP diffusr_heat_diffusion_(SEXP v0SEXP, SEXP WSEXP, SEXP bSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(diffusr_heat_diffusion__try(v0SEXP, WSEXP, bSEXP));
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
static SEXP diffusr_stoch_col_norm__try(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(stoch_col_norm_(W));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP diffusr_stoch_col_norm_(SEXP WSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(diffusr_stoch_col_norm__try(WSEXP));
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
static SEXP diffusr_laplacian__try(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(laplacian_(W));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP diffusr_laplacian_(SEXP WSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(diffusr_laplacian__try(WSEXP));
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
Eigen::VectorXd mrwr_(const Eigen::VectorXd& p0, const Eigen::MatrixXd& W, const double r);
static SEXP diffusr_mrwr__try(SEXP p0SEXP, SEXP WSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(mrwr_(p0, W, r));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP diffusr_mrwr_(SEXP p0SEXP, SEXP WSEXP, SEXP rSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(diffusr_mrwr__try(p0SEXP, WSEXP, rSEXP));
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
Rcpp::List neighbors_(const Rcpp::IntegerVector& node_idxs, const Rcpp::NumericMatrix& W, const int k, const bool use_edge_weights);
static SEXP diffusr_neighbors__try(SEXP node_idxsSEXP, SEXP WSEXP, SEXP kSEXP, SEXP use_edge_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type node_idxs(node_idxsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_edge_weights(use_edge_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(neighbors_(node_idxs, W, k, use_edge_weights));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP diffusr_neighbors_(SEXP node_idxsSEXP, SEXP WSEXP, SEXP kSEXP, SEXP use_edge_weightsSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(diffusr_neighbors__try(node_idxsSEXP, WSEXP, kSEXP, use_edge_weightsSEXP));
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
static int diffusr_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Eigen::VectorXd(*.heat_diffusion.cpp)(const Eigen::VectorXd&,const Eigen::MatrixXd&,const double)");
        signatures.insert("Eigen::MatrixXd(*.stoch.col.norm.cpp)(const Eigen::MatrixXd&)");
        signatures.insert("Eigen::MatrixXd(*.laplacian.cpp)(const Eigen::MatrixXd&)");
        signatures.insert("Eigen::VectorXd(*.mrwr.cpp)(const Eigen::VectorXd&,const Eigen::MatrixXd&,const double)");
        signatures.insert("Rcpp::List(*.neighbors.cpp)(const Rcpp::IntegerVector&,const Rcpp::NumericMatrix&,const int,const bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP diffusr_RcppExport_registerCCallable() { 
    R_RegisterCCallable("diffusr", "diffusr_.heat_diffusion.cpp", (DL_FUNC)diffusr_heat_diffusion__try);
    R_RegisterCCallable("diffusr", "diffusr_.stoch.col.norm.cpp", (DL_FUNC)diffusr_stoch_col_norm__try);
    R_RegisterCCallable("diffusr", "diffusr_.laplacian.cpp", (DL_FUNC)diffusr_laplacian__try);
    R_RegisterCCallable("diffusr", "diffusr_.mrwr.cpp", (DL_FUNC)diffusr_mrwr__try);
    R_RegisterCCallable("diffusr", "diffusr_.neighbors.cpp", (DL_FUNC)diffusr_neighbors__try);
    R_RegisterCCallable("diffusr", "diffusr_RcppExport_validate", (DL_FUNC)diffusr_RcppExport_validate);
    return R_NilValue;
}
