// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ecor
Eigen::MatrixXd ecor(Eigen::MatrixXd mat);
RcppExport SEXP _MeDuSAJ_ecor(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(ecor(mat));
    return rcpp_result_gen;
END_RCPP
}
// reml
std::vector<Eigen::MatrixXd> reml(Eigen::VectorXd start, Eigen::MatrixXd& X, Eigen::VectorXd& y, std::vector<Eigen::MatrixXd>& Z, int maxiter);
RcppExport SEXP _MeDuSAJ_reml(SEXP startSEXP, SEXP XSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type start(startSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<Eigen::MatrixXd>& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(reml(start, X, y, Z, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// reml2
std::vector<Eigen::MatrixXd> reml2(Eigen::VectorXd start, Eigen::MatrixXd& X, Eigen::VectorXd& y, std::vector<Eigen::MatrixXd>& Z, int maxiter);
RcppExport SEXP _MeDuSAJ_reml2(SEXP startSEXP, SEXP XSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type start(startSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::vector<Eigen::MatrixXd>& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(reml2(start, X, y, Z, maxiter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MeDuSAJ_ecor", (DL_FUNC) &_MeDuSAJ_ecor, 1},
    {"_MeDuSAJ_reml", (DL_FUNC) &_MeDuSAJ_reml, 5},
    {"_MeDuSAJ_reml2", (DL_FUNC) &_MeDuSAJ_reml2, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MeDuSAJ(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
