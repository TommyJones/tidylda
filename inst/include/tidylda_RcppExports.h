// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_tidylda_RCPPEXPORTS_H_GEN_
#define RCPP_tidylda_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace tidylda {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("tidylda", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("tidylda", "_tidylda_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in tidylda");
            }
        }
    }

    inline Rcpp::List create_lexicon(const IntegerMatrix& Cd_in, const NumericMatrix& Beta_in, const arma::sp_mat& dtm_in, const std::vector<double>& alpha, const bool& freeze_topics) {
        typedef SEXP(*Ptr_create_lexicon)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_create_lexicon p_create_lexicon = NULL;
        if (p_create_lexicon == NULL) {
            validateSignature("Rcpp::List(*create_lexicon)(const IntegerMatrix&,const NumericMatrix&,const arma::sp_mat&,const std::vector<double>&,const bool&)");
            p_create_lexicon = (Ptr_create_lexicon)R_GetCCallable("tidylda", "_tidylda_create_lexicon");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_create_lexicon(Shield<SEXP>(Rcpp::wrap(Cd_in)), Shield<SEXP>(Rcpp::wrap(Beta_in)), Shield<SEXP>(Rcpp::wrap(dtm_in)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(freeze_topics)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List fit_lda_c(const std::vector<std::vector<std::size_t>>& Docs, const std::vector<std::vector<std::size_t>>& Zd_in, const IntegerMatrix& Cd_in, const IntegerMatrix& Cv_in, const std::vector<long>& Ck_in, const std::vector<double> alpha_in, const NumericMatrix& eta_in, const std::size_t& iterations, const int& burnin, const bool& optimize_alpha, const bool& calc_likelihood, const NumericMatrix& Beta_in, const bool& freeze_topics, const std::size_t& threads = 1, const bool& verbose = false) {
        typedef SEXP(*Ptr_fit_lda_c)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_fit_lda_c p_fit_lda_c = NULL;
        if (p_fit_lda_c == NULL) {
            validateSignature("Rcpp::List(*fit_lda_c)(const std::vector<std::vector<std::size_t>>&,const std::vector<std::vector<std::size_t>>&,const IntegerMatrix&,const IntegerMatrix&,const std::vector<long>&,const std::vector<double>,const NumericMatrix&,const std::size_t&,const int&,const bool&,const bool&,const NumericMatrix&,const bool&,const std::size_t&,const bool&)");
            p_fit_lda_c = (Ptr_fit_lda_c)R_GetCCallable("tidylda", "_tidylda_fit_lda_c");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_fit_lda_c(Shield<SEXP>(Rcpp::wrap(Docs)), Shield<SEXP>(Rcpp::wrap(Zd_in)), Shield<SEXP>(Rcpp::wrap(Cd_in)), Shield<SEXP>(Rcpp::wrap(Cv_in)), Shield<SEXP>(Rcpp::wrap(Ck_in)), Shield<SEXP>(Rcpp::wrap(alpha_in)), Shield<SEXP>(Rcpp::wrap(eta_in)), Shield<SEXP>(Rcpp::wrap(iterations)), Shield<SEXP>(Rcpp::wrap(burnin)), Shield<SEXP>(Rcpp::wrap(optimize_alpha)), Shield<SEXP>(Rcpp::wrap(calc_likelihood)), Shield<SEXP>(Rcpp::wrap(Beta_in)), Shield<SEXP>(Rcpp::wrap(freeze_topics)), Shield<SEXP>(Rcpp::wrap(threads)), Shield<SEXP>(Rcpp::wrap(verbose)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

}

#endif // RCPP_tidylda_RCPPEXPORTS_H_GEN_
