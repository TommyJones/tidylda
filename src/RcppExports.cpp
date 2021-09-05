// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// create_lexicon
Rcpp::List create_lexicon(const IntegerMatrix& Cd_in, const NumericMatrix& Beta_in, const arma::sp_mat& dtm_in, const std::vector<double>& alpha, const bool& freeze_topics);
static SEXP _tidylda_create_lexicon_try(SEXP Cd_inSEXP, SEXP Beta_inSEXP, SEXP dtm_inSEXP, SEXP alphaSEXP, SEXP freeze_topicsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type Cd_in(Cd_inSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Beta_in(Beta_inSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type dtm_in(dtm_inSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type freeze_topics(freeze_topicsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_lexicon(Cd_in, Beta_in, dtm_in, alpha, freeze_topics));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidylda_create_lexicon(SEXP Cd_inSEXP, SEXP Beta_inSEXP, SEXP dtm_inSEXP, SEXP alphaSEXP, SEXP freeze_topicsSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidylda_create_lexicon_try(Cd_inSEXP, Beta_inSEXP, dtm_inSEXP, alphaSEXP, freeze_topicsSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
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
// fit_lda_c
Rcpp::List fit_lda_c(const std::vector<std::vector<std::size_t>>& Docs, const std::vector<std::vector<std::size_t>>& Zd_in, const IntegerMatrix& Cd_in, const IntegerMatrix& Cv_in, const std::vector<long>& Ck_in, const std::vector<double> alpha_in, const NumericMatrix& eta_in, const std::size_t& iterations, const int& burnin, const bool& optimize_alpha, const bool& calc_likelihood, const NumericMatrix& Beta_in, const bool& freeze_topics, const std::size_t& threads, const bool& verbose);
static SEXP _tidylda_fit_lda_c_try(SEXP DocsSEXP, SEXP Zd_inSEXP, SEXP Cd_inSEXP, SEXP Cv_inSEXP, SEXP Ck_inSEXP, SEXP alpha_inSEXP, SEXP eta_inSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP optimize_alphaSEXP, SEXP calc_likelihoodSEXP, SEXP Beta_inSEXP, SEXP freeze_topicsSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::size_t>>& >::type Docs(DocsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<std::size_t>>& >::type Zd_in(Zd_inSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type Cd_in(Cd_inSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type Cv_in(Cv_inSEXP);
    Rcpp::traits::input_parameter< const std::vector<long>& >::type Ck_in(Ck_inSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type alpha_in(alpha_inSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type eta_in(eta_inSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< const int& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool& >::type optimize_alpha(optimize_alphaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type calc_likelihood(calc_likelihoodSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Beta_in(Beta_inSEXP);
    Rcpp::traits::input_parameter< const bool& >::type freeze_topics(freeze_topicsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_lda_c(Docs, Zd_in, Cd_in, Cv_in, Ck_in, alpha_in, eta_in, iterations, burnin, optimize_alpha, calc_likelihood, Beta_in, freeze_topics, threads, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _tidylda_fit_lda_c(SEXP DocsSEXP, SEXP Zd_inSEXP, SEXP Cd_inSEXP, SEXP Cv_inSEXP, SEXP Ck_inSEXP, SEXP alpha_inSEXP, SEXP eta_inSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP optimize_alphaSEXP, SEXP calc_likelihoodSEXP, SEXP Beta_inSEXP, SEXP freeze_topicsSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_tidylda_fit_lda_c_try(DocsSEXP, Zd_inSEXP, Cd_inSEXP, Cv_inSEXP, Ck_inSEXP, alpha_inSEXP, eta_inSEXP, iterationsSEXP, burninSEXP, optimize_alphaSEXP, calc_likelihoodSEXP, Beta_inSEXP, freeze_topicsSEXP, threadsSEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
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
static int _tidylda_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*create_lexicon)(const IntegerMatrix&,const NumericMatrix&,const arma::sp_mat&,const std::vector<double>&,const bool&)");
        signatures.insert("Rcpp::List(*fit_lda_c)(const std::vector<std::vector<std::size_t>>&,const std::vector<std::vector<std::size_t>>&,const IntegerMatrix&,const IntegerMatrix&,const std::vector<long>&,const std::vector<double>,const NumericMatrix&,const std::size_t&,const int&,const bool&,const bool&,const NumericMatrix&,const bool&,const std::size_t&,const bool&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _tidylda_RcppExport_registerCCallable() { 
    R_RegisterCCallable("tidylda", "_tidylda_create_lexicon", (DL_FUNC)_tidylda_create_lexicon_try);
    R_RegisterCCallable("tidylda", "_tidylda_fit_lda_c", (DL_FUNC)_tidylda_fit_lda_c_try);
    R_RegisterCCallable("tidylda", "_tidylda_RcppExport_validate", (DL_FUNC)_tidylda_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_tidylda_create_lexicon", (DL_FUNC) &_tidylda_create_lexicon, 5},
    {"_tidylda_fit_lda_c", (DL_FUNC) &_tidylda_fit_lda_c, 15},
    {"_tidylda_RcppExport_registerCCallable", (DL_FUNC) &_tidylda_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_tidylda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
