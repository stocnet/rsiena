// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

// Not needed currently, because flattening function for one chain now in siena07Utilities
// // [[Rcpp::export]]
// DataFrame flattenContribList(List chainsList) {
//   std::vector<int> chain_col, group_col, period_col, ministep_col, effect_col, choice_col;
//   std::vector<std::string> effectname_col, effecttype_col;
//   std::vector<double> value_col;

//   int nChains = chainsList.size();
//   for (int ch = 0; ch < nChains; ++ch) {
//     List groupList = chainsList[ch];       // [group]
//     int nGroups = groupList.size();
//     for (int g = 0; g < nGroups; ++g) {
//       List periodList = groupList[g];      // [period]
//       int nPeriods = periodList.size();
//       for (int p = 0; p < nPeriods; ++p) {
//         List ministeps = periodList[p];
//         int nMinisteps = ministeps.size();
//         for (int m = 0; m < nMinisteps; ++m) {
//           RObject mat_obj = ministeps[m];
//           if (mat_obj.isNULL() || !Rf_isMatrix(mat_obj) || TYPEOF(mat_obj) != REALSXP) continue;
//           NumericMatrix mat(mat_obj);
//           CharacterVector effectNames = mat.attr("effectNames");
//           CharacterVector effectTypes = mat.attr("effectTypes");
//           int nEff = mat.nrow(), nChoice = mat.ncol();
//           for (int e = 0; e < nEff; ++e) {
//             std::string effname = (effectNames.size() > 0) ? as<std::string>(effectNames[e]) : "";
//             std::string efftype = (effectTypes.size() > 0) ? as<std::string>(effectTypes[e]) : "";
//             for (int c = 0; c < nChoice; ++c) {
//               chain_col.push_back(ch+1);
//               group_col.push_back(g+1);
//               period_col.push_back(p+1);
//               ministep_col.push_back(m+1);
//               effect_col.push_back(e+1);
//               choice_col.push_back(c+1);
//               effectname_col.push_back(effname);
//               effecttype_col.push_back(efftype);
//               value_col.push_back(mat(e,c));
//             }
//           }
//         }
//       }
//     }
//   }
//   return DataFrame::create(
//     _["chain"]      = chain_col,
//     _["group"]      = group_col,
//     _["period"]     = period_col,
//     _["ministep"]   = ministep_col,
//     _["effect"]     = effect_col,
//     _["choice"]     = choice_col,
//     _["effectname"] = effectname_col,
//     _["effecttype"] = effecttype_col,
//     _["contribution"]= value_col,
//     _["stringsAsFactors"]= false
//   );
// }

static void softmax_inplace(const double* x, double* out, int n) {
  double mx = *std::max_element(x, x + n);
  double s = 0.0;
  for (int i = 0; i < n; i++) { out[i] = std::exp(x[i] - mx); s += out[i]; }
  for (int i = 0; i < n; i++) out[i] /= s;
}

// [[Rcpp::export]]
arma::vec softmax_arma(const arma::vec& x) {
  arma::vec out(x.n_elem);
  softmax_inplace(x.memptr(), out.memptr(), x.n_elem);
  return out;
}

// [[Rcpp::export]]
arma::vec softmax_arma_by_group(const arma::vec& x, const arma::ivec& group) {
  int n = x.n_elem;
  arma::vec out(n);          // one allocation for everything
  int start = 0;
  while (start < n) {
    int g = group[start], end = start + 1;
    while (end < n && group[end] == g) end++;           // find group boundary
    softmax_inplace(x.memptr() + start,                 // pointer into input slice
                    out.memptr() + start,               // pointer into output slice
                    end - start);                       // group size
    start = end;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec softmax_rcpp_grouped(const arma::vec& x,
                                const arma::ivec& g1,
                                const arma::ivec& g2,
                                const arma::ivec& g3) {
  int n = x.n_elem;
  arma::vec out(n);
  int start = 0;
  while (start < n) {
    int end = start + 1;
    while (end < n && g1[end] == g1[start] &&
                      g2[end] == g2[start] &&
                      g3[end] == g3[start]) end++;
    softmax_inplace(x.memptr() + start, out.memptr() + start, end - start);
    start = end;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec softmax_rcpp_grouped_mat(const arma::vec& x, const arma::imat& G) {
  int n = x.n_elem, ncols = G.n_cols;
  arma::vec out(n);
  std::vector<const int*> cols(ncols);
  for (int c = 0; c < ncols; c++) cols[c] = G.colptr(c);  // contiguous per column
  int start = 0;
  while (start < n) {
    int end = start + 1;
    while (end < n) {
      bool same = true;
      for (int c = 0; c < ncols && same; c++)
        if (cols[c][end] != cols[c][start]) same = false;
      if (!same) break;
      end++;
    }
    softmax_inplace(x.memptr() + start, out.memptr() + start, end - start);
    start = end;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec softmax_rcpp_grouped_lst(const arma::vec& x, const Rcpp::List& G) {
  int n = x.n_elem, ncols = G.size();
  arma::vec out(n);
  std::vector<const int*> cols(ncols);
  for (int c = 0; c < ncols; c++) {
    Rcpp::IntegerVector v = G[c];   // zero-copy wrap of existing R SEXP data
    cols[c] = v.begin();
  }
  int start = 0;
  while (start < n) {
    int end = start + 1;
    while (end < n) {
      bool same = true;
      for (int c = 0; c < ncols && same; c++)
        if (cols[c][end] != cols[c][start]) same = false;
      if (!same) break;
      end++;
    }
    softmax_inplace(x.memptr() + start, out.memptr() + start, end - start);
    start = end;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec softmax_rcpp_grouped_cols(const Rcpp::DataFrame& data,
                                   const std::string& val_col,
                                   const Rcpp::CharacterVector& group_cols) {
  Rcpp::NumericVector x_r = data[val_col];
  const double* x = x_r.begin();
  int n = x_r.size(), ncols = group_cols.size();
  arma::vec out(n);
  std::vector<const int*> cols(ncols);
  for (int c = 0; c < ncols; c++) {
    Rcpp::IntegerVector v = data[Rcpp::as<std::string>(group_cols[c])];
    cols[c] = v.begin();
  }
  int start = 0;
  while (start < n) {
    int end = start + 1;
    while (end < n) {
      bool same = true;
      for (int c = 0; c < ncols && same; c++)
        if (cols[c][end] != cols[c][start]) same = false;
      if (!same) break;
      end++;
    }
    softmax_inplace(x + start, out.memptr() + start, end - start);
    start = end;
  }
  return out;
}