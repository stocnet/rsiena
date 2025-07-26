// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>

using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
arma::vec softmax_arma(const arma::vec& x) {
  double max_x = x.max();
  arma::vec expx = arma::exp(x - max_x);
  return expx / arma::sum(expx);
}

// [[Rcpp::export]]
arma::vec softmax_arma_by_group(const arma::vec& x, const arma::ivec& group) {
  int n = x.n_elem;
  arma::vec out(n);
  std::unordered_map<int, std::vector<int>> group_map;
  for (int i = 0; i < n; ++i)
    group_map[group[i]].push_back(i);

  for (auto& kv : group_map) {
    auto& idxs = kv.second;
    arma::vec subv(idxs.size());
    for (size_t j = 0; j < idxs.size(); ++j)
      subv[j] = x[idxs[j]];
    double m = subv.max();
    arma::vec exps = arma::exp(subv - m);
    double sum_exps = arma::sum(exps);
    for (size_t j = 0; j < idxs.size(); ++j)
      out[idxs[j]] = exps[j] / sum_exps;
  }

  return out;
}
