// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4exposure_model_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4outcome_model_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4outcome_model_continuous_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_restrictBeta_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4outcome_model_restrictBeta_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4exposure_model_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4exposure_model_mod, 0},
    {"_rcpp_module_boot_stan_fit4outcome_model_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4outcome_model_mod, 0},
    {"_rcpp_module_boot_stan_fit4outcome_model_continuous_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4outcome_model_continuous_mod, 0},
    {"_rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_mod, 0},
    {"_rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_restrictBeta_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4outcome_model_multipleExpCurves_restrictBeta_mod, 0},
    {"_rcpp_module_boot_stan_fit4outcome_model_restrictBeta_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4outcome_model_restrictBeta_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bercs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}