#include <Rcpp.h>
using namespace Rcpp;

/**
 * This set of functions implements division with a single modification: 0/0 --> 0.
 * This might not be mathematically "correct", but occurrs when calculating the
 * residuals for genes where all y = 0.
 *
 * It is roughly half as fast as the native `/` operator, but this is still much better
 * than having to create a mask to check for those cases afterwards.
 */

// [[Rcpp::export(rng = false)]]
NumericVector div_zbz_dbl(const NumericVector a, const NumericVector b) {
  size_t as = a.size();
  size_t bs = b.size();
  if (as != bs) {
    stop("Size of a and b must match");
  }
  NumericVector res(as);
  for (size_t idx = 0; idx < as; idx++) {
    double ai = a[idx];
    double bi = b[idx];
    if (ai == 0 && bi == 0) {
      res[idx] = 0;
    } else {
      res[idx] = ai / bi;
    }
  }
  return res;
}

// [[Rcpp::export(rng = false)]]
NumericVector div_zbz_int(const IntegerVector a, const IntegerVector b) {
  size_t as = a.size();
  size_t bs = b.size();
  if (as != bs) {
    stop("Size of a and b must match");
  }
  NumericVector res(as);
  for (size_t idx = 0; idx < as; idx++) {
    size_t ai = a[idx];
    size_t bi = b[idx];
    if (ai == 0 && bi == 0) {
      res[idx] = 0;
    } else {
      res[idx] = (double)ai / bi;
    }
  }
  return res;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix div_zbz_dbl_mat(const NumericMatrix a, const NumericMatrix b) {
  if (a.nrow() != b.nrow() || a.ncol() != b.ncol()) {
    stop("The dimensions of the matrices must match");
  }
  NumericVector vec = div_zbz_dbl(a, b);
  NumericMatrix res(a.nrow(), a.ncol(), vec.begin());
  return res;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix div_zbz_int_mat(const IntegerMatrix a, const IntegerMatrix b) {
  if (a.nrow() != b.nrow() || a.ncol() != b.ncol()) {
    stop("The dimensions of the matrices must match");
  }
  NumericVector vec = div_zbz_int(a, b);
  NumericMatrix res(a.nrow(), a.ncol(), vec.begin());
  return res;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
div_zbz(0:4, 0:4)
*/
