// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_diffusr_RCPPEXPORTS_H_GEN_
#define RCPP_diffusr_RCPPEXPORTS_H_GEN_

#include <RcppEigen.h>
#include <Rcpp.h>

namespace diffusr {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("diffusr", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("diffusr", "_diffusr_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in diffusr");
            }
        }
    }

    inline MatrixXd heat_diffusion_(const MatrixXd& v0, const MatrixXd& W, const double t) {
        typedef SEXP(*Ptr_heat_diffusion_)(SEXP,SEXP,SEXP);
        static Ptr_heat_diffusion_ p_heat_diffusion_ = NULL;
        if (p_heat_diffusion_ == NULL) {
            validateSignature("MatrixXd(*heat_diffusion_)(const MatrixXd&,const MatrixXd&,const double)");
            p_heat_diffusion_ = (Ptr_heat_diffusion_)R_GetCCallable("diffusr", "_diffusr_heat_diffusion_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_heat_diffusion_(Shield<SEXP>(Rcpp::wrap(v0)), Shield<SEXP>(Rcpp::wrap(W)), Shield<SEXP>(Rcpp::wrap(t)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline MatrixXd stoch_col_norm_(const MatrixXd& W) {
        typedef SEXP(*Ptr_stoch_col_norm_)(SEXP);
        static Ptr_stoch_col_norm_ p_stoch_col_norm_ = NULL;
        if (p_stoch_col_norm_ == NULL) {
            validateSignature("MatrixXd(*stoch_col_norm_)(const MatrixXd&)");
            p_stoch_col_norm_ = (Ptr_stoch_col_norm_)R_GetCCallable("diffusr", "_diffusr_stoch_col_norm_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_stoch_col_norm_(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline SpMat stoch_col_norm_s(const SpMat& W) {
        typedef SEXP(*Ptr_stoch_col_norm_s)(SEXP);
        static Ptr_stoch_col_norm_s p_stoch_col_norm_s = NULL;
        if (p_stoch_col_norm_s == NULL) {
            validateSignature("SpMat(*stoch_col_norm_s)(const SpMat&)");
            p_stoch_col_norm_s = (Ptr_stoch_col_norm_s)R_GetCCallable("diffusr", "_diffusr_stoch_col_norm_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_stoch_col_norm_s(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SpMat >(rcpp_result_gen);
    }

    inline MatrixXd laplacian_(const MatrixXd& W) {
        typedef SEXP(*Ptr_laplacian_)(SEXP);
        static Ptr_laplacian_ p_laplacian_ = NULL;
        if (p_laplacian_ == NULL) {
            validateSignature("MatrixXd(*laplacian_)(const MatrixXd&)");
            p_laplacian_ = (Ptr_laplacian_)R_GetCCallable("diffusr", "_diffusr_laplacian_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_laplacian_(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline MatrixXd laplacian_s(const SpMat& W) {
        typedef SEXP(*Ptr_laplacian_s)(SEXP);
        static Ptr_laplacian_s p_laplacian_s = NULL;
        if (p_laplacian_s == NULL) {
            validateSignature("MatrixXd(*laplacian_s)(const SpMat&)");
            p_laplacian_s = (Ptr_laplacian_s)R_GetCCallable("diffusr", "_diffusr_laplacian_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_laplacian_s(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline VectorXd node_degrees_(const MatrixXd& W) {
        typedef SEXP(*Ptr_node_degrees_)(SEXP);
        static Ptr_node_degrees_ p_node_degrees_ = NULL;
        if (p_node_degrees_ == NULL) {
            validateSignature("VectorXd(*node_degrees_)(const MatrixXd&)");
            p_node_degrees_ = (Ptr_node_degrees_)R_GetCCallable("diffusr", "_diffusr_node_degrees_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_node_degrees_(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<VectorXd >(rcpp_result_gen);
    }

    inline VectorXd node_degrees_s(const MSpMat& W) {
        typedef SEXP(*Ptr_node_degrees_s)(SEXP);
        static Ptr_node_degrees_s p_node_degrees_s = NULL;
        if (p_node_degrees_s == NULL) {
            validateSignature("VectorXd(*node_degrees_s)(const MSpMat&)");
            p_node_degrees_s = (Ptr_node_degrees_s)R_GetCCallable("diffusr", "_diffusr_node_degrees_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_node_degrees_s(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<VectorXd >(rcpp_result_gen);
    }

    inline MatrixXd hub_normalize_(const MatrixXd& W) {
        typedef SEXP(*Ptr_hub_normalize_)(SEXP);
        static Ptr_hub_normalize_ p_hub_normalize_ = NULL;
        if (p_hub_normalize_ == NULL) {
            validateSignature("MatrixXd(*hub_normalize_)(const MatrixXd&)");
            p_hub_normalize_ = (Ptr_hub_normalize_)R_GetCCallable("diffusr", "_diffusr_hub_normalize_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hub_normalize_(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline MatrixXd hub_normalize_s(const SpMat& W) {
        typedef SEXP(*Ptr_hub_normalize_s)(SEXP);
        static Ptr_hub_normalize_s p_hub_normalize_s = NULL;
        if (p_hub_normalize_s == NULL) {
            validateSignature("MatrixXd(*hub_normalize_s)(const SpMat&)");
            p_hub_normalize_s = (Ptr_hub_normalize_s)R_GetCCallable("diffusr", "_diffusr_hub_normalize_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hub_normalize_s(Shield<SEXP>(Rcpp::wrap(W)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<MatrixXd >(rcpp_result_gen);
    }

    inline VectorXd mrwr_(const MatrixXd& p0, const MatrixXd& W, const double r, const double thresh, const int niter, const bool do_analytical) {
        typedef SEXP(*Ptr_mrwr_)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mrwr_ p_mrwr_ = NULL;
        if (p_mrwr_ == NULL) {
            validateSignature("VectorXd(*mrwr_)(const MatrixXd&,const MatrixXd&,const double,const double,const int,const bool)");
            p_mrwr_ = (Ptr_mrwr_)R_GetCCallable("diffusr", "_diffusr_mrwr_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mrwr_(Shield<SEXP>(Rcpp::wrap(p0)), Shield<SEXP>(Rcpp::wrap(W)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(niter)), Shield<SEXP>(Rcpp::wrap(do_analytical)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<VectorXd >(rcpp_result_gen);
    }

    inline VectorXd mrwr_s(const MatrixXd& p0, const SpMat& W, const double r, const double thresh, const int niter, const bool do_analytical) {
        typedef SEXP(*Ptr_mrwr_s)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mrwr_s p_mrwr_s = NULL;
        if (p_mrwr_s == NULL) {
            validateSignature("VectorXd(*mrwr_s)(const MatrixXd&,const SpMat&,const double,const double,const int,const bool)");
            p_mrwr_s = (Ptr_mrwr_s)R_GetCCallable("diffusr", "_diffusr_mrwr_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mrwr_s(Shield<SEXP>(Rcpp::wrap(p0)), Shield<SEXP>(Rcpp::wrap(W)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(thresh)), Shield<SEXP>(Rcpp::wrap(niter)), Shield<SEXP>(Rcpp::wrap(do_analytical)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<VectorXd >(rcpp_result_gen);
    }

    inline List neighbors_(const vector<int>& node_idxs, const MatrixXd& W, const int& k) {
        typedef SEXP(*Ptr_neighbors_)(SEXP,SEXP,SEXP);
        static Ptr_neighbors_ p_neighbors_ = NULL;
        if (p_neighbors_ == NULL) {
            validateSignature("List(*neighbors_)(const vector<int>&,const MatrixXd&,const int&)");
            p_neighbors_ = (Ptr_neighbors_)R_GetCCallable("diffusr", "_diffusr_neighbors_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_neighbors_(Shield<SEXP>(Rcpp::wrap(node_idxs)), Shield<SEXP>(Rcpp::wrap(W)), Shield<SEXP>(Rcpp::wrap(k)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List neighbors_s(const vector<int>& node_idxs, const MSpMat& W, const int& k) {
        typedef SEXP(*Ptr_neighbors_s)(SEXP,SEXP,SEXP);
        static Ptr_neighbors_s p_neighbors_s = NULL;
        if (p_neighbors_s == NULL) {
            validateSignature("List(*neighbors_s)(const vector<int>&,const MSpMat&,const int&)");
            p_neighbors_s = (Ptr_neighbors_s)R_GetCCallable("diffusr", "_diffusr_neighbors_s");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_neighbors_s(Shield<SEXP>(Rcpp::wrap(node_idxs)), Shield<SEXP>(Rcpp::wrap(W)), Shield<SEXP>(Rcpp::wrap(k)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_diffusr_RCPPEXPORTS_H_GEN_
