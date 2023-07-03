#include <Rcpp.h>
#include "Rtatami.h"
#include <vector>
#include <algorithm>

using namespace Rcpp;


// [[Rcpp::export(rng=false)]]
NumericVector column_sums(RObject initmat){
  Rtatami::BoundNumericPointer parsed(initmat);

  const auto& ptr = parsed->ptr;

  auto NR = ptr->nrow();
  auto NC = ptr->ncol();
  std::vector<double> buffer(NR);
  NumericVector output(NC);
  auto wrk = ptr->dense_column();

  for(int i = 0; i < NC; ++i) {
    auto extracted = wrk->fetch(i, buffer.data());
    output[i] = std::accumulate(extracted, extracted + NR, 0.0);
  }

  return output;
}




