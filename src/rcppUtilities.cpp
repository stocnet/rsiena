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

// // [[Rcpp::export]]
// Rcpp::List flattenChangeContributionsWide_rcpp(
//     Rcpp::List chains,
//     Rcpp::CharacterVector effectNames,
//     Rcpp::Nullable<Rcpp::CharacterVector> depvar_sxp = R_NilValue)
// {
//     int nEffSel = effectNames.size();
//     int nChains = chains.size();

//     bool filterDepvar = depvar_sxp.isNotNull();
//     Rcpp::CharacterVector depvar;
//     if (filterDepvar) depvar = depvar_sxp;

//     auto keepNet = [&](Rcpp::RObject obj) -> bool {
//         if (!filterDepvar) return true;
//         if (!obj.hasAttribute("networkName")) return false;
//         std::string net = Rcpp::as<std::string>(obj.attr("networkName"));
//         for (int d = 0; d < (int)depvar.size(); d++)
//             if (net == Rcpp::as<std::string>(depvar[d])) return true;
//         return false;
//     };

//     /* discover raw effect order from first non-null kept ministep */
//     Rcpp::CharacterVector rawEffNames;
//     bool foundRaw = false;
//     for (int ch = 0; ch < nChains && !foundRaw; ch++) {
//         Rcpp::List cl = chains[ch];
//         for (int g = 0; g < cl.size() && !foundRaw; g++) {
//             Rcpp::List gl = cl[g];
//             for (int p = 0; p < gl.size() && !foundRaw; p++) {
//                 Rcpp::List pl = gl[p];
//                 for (int m = 0; m < pl.size(); m++) {
//                     Rcpp::RObject obj = pl[m];
//                     if (obj.isNULL() || !keepNet(obj)) continue;
//                     if (obj.hasAttribute("effectNames")) {
//                         rawEffNames = obj.attr("effectNames");
//                         foundRaw = true;
//                     }
//                     break;
//                 }
//             }
//         }
//     }

//     /* matchIdx[e] = row in raw mat for effectNames[e] */
//     std::vector<int> matchIdx(nEffSel, -1);
//     for (int e = 0; e < nEffSel; e++) {
//         std::string want = Rcpp::as<std::string>(effectNames[e]);
//         for (int r = 0; r < (int)rawEffNames.size(); r++)
//             if (want == Rcpp::as<std::string>(rawEffNames[r]))
//                 { matchIdx[e] = r; break; }
//     }

//     /* count pass */
//     int totalRows = 0;
//     for (int ch = 0; ch < nChains; ch++) {
//         Rcpp::List cl = chains[ch];
//         for (int g = 0; g < cl.size(); g++) {
//             Rcpp::List gl = cl[g];
//             for (int p = 0; p < gl.size(); p++) {
//                 Rcpp::List pl = gl[p];
//                 for (int m = 0; m < pl.size(); m++) {
//                     Rcpp::RObject obj = pl[m];
//                     if (obj.isNULL() || !keepNet(obj)) continue;
//                     totalRows += Rcpp::NumericMatrix(obj).ncol();
//                 }
//             }
//         }
//     }

//     /* allocate */
//     Rcpp::NumericMatrix contrib(totalRows, nEffSel);
//     Rcpp::IntegerVector chain_v(totalRows), group_v(totalRows),
//                         period_v(totalRows), mstep_v(totalRows),
//                         choice_v(totalRows), grpid_v(totalRows);

//     /* fill pass */
//     int row = 0, grp = 0;
//     for (int ch = 0; ch < nChains; ch++) {
//         Rcpp::List cl = chains[ch];
//         for (int g = 0; g < cl.size(); g++) {
//             Rcpp::List gl = cl[g];
//             for (int p = 0; p < gl.size(); p++) {
//                 Rcpp::List pl = gl[p];
//                 for (int m = 0; m < pl.size(); m++) {
//                     Rcpp::RObject obj = pl[m];
//                     if (obj.isNULL() || !keepNet(obj)) continue;
//                     Rcpp::NumericMatrix mat(obj);
//                     int nEff = mat.nrow(), nChoice = mat.ncol();
//                     grp++;
//                     for (int c = 0; c < nChoice; c++) {
//                         for (int e = 0; e < nEffSel; e++)
//                           contrib(row, e) = (matchIdx[e] >= 0 && matchIdx[e] < nEff)
//                               ? mat(matchIdx[e], c) : NA_REAL;
//                         chain_v[row] = ch + 1;  group_v[row] = g + 1;
//                         period_v[row] = p + 1;  mstep_v[row] = m + 1;
//                         choice_v[row] = c + 1;  grpid_v[row] = grp;
//                         row++;
//                     }
//                 }
//             }
//         }
//     }

//     Rcpp::colnames(contrib) = effectNames;

//     return Rcpp::List::create(
//         Rcpp::Named("contribMat") = contrib,
//         Rcpp::Named("chain")       = chain_v,
//         Rcpp::Named("group")       = group_v,
//         Rcpp::Named("period")      = period_v,
//         Rcpp::Named("ministep")    = mstep_v,
//         Rcpp::Named("choice")      = choice_v,
//         Rcpp::Named("group_id")    = grpid_v
//     );
// }

// ---------------------------------------------------------------------------
// Grouped aggregation: single-pass mean (or sum) over contiguous groups
// defined by an integer matrix of group keys.
// x:       numeric vector of values
// G:       integer matrix (nrow = length(x), ncol = number of group columns)
//          Rows must be sorted so that identical key rows are contiguous.
// na_rm:   if true, NA values are skipped
// do_mean: if true, compute mean; if false, compute sum
// Returns a List with:
//   $value   numeric vector of group means/sums
//   $count   integer vector of group sizes (excluding NAs if na_rm)
//   $key     integer matrix of unique group keys (nrow = nGroups)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grouped_agg_cpp(const Rcpp::NumericVector& x,
                           const Rcpp::IntegerMatrix& G,
                           bool na_rm = true,
                           bool do_mean = true) {
    int n = x.size(), ncols = G.ncol();
    if (G.nrow() != n)
        Rcpp::stop("length(x) must equal nrow(G)");

    // Pointers into each column of G for fast comparison
    std::vector<const int*> cols(ncols);
    for (int c = 0; c < ncols; c++) cols[c] = &G(0, c);

    // Count groups (single pass)
    int nG = 0;
    {
        int start = 0;
        while (start < n) {
            nG++;
            int end = start + 1;
            while (end < n) {
                bool same = true;
                for (int c = 0; c < ncols && same; c++)
                    if (cols[c][end] != cols[c][start]) same = false;
                if (!same) break;
                end++;
            }
            start = end;
        }
    }

    // Allocate output
    Rcpp::NumericVector val(nG);
    Rcpp::IntegerVector cnt(nG);
    Rcpp::IntegerMatrix key(nG, ncols);

    // Fill pass
    int gi = 0, start = 0;
    while (start < n) {
        // Store group key
        for (int c = 0; c < ncols; c++) key(gi, c) = cols[c][start];

        // Find group end
        int end = start + 1;
        while (end < n) {
            bool same = true;
            for (int c = 0; c < ncols && same; c++)
                if (cols[c][end] != cols[c][start]) same = false;
            if (!same) break;
            end++;
        }

        // Accumulate
        double sum = 0.0;
        int count = 0;
        for (int i = start; i < end; i++) {
            if (na_rm && ISNAN(x[i])) continue;
            sum += x[i];
            count++;
        }
        val[gi] = do_mean ? (count > 0 ? sum / count : NA_REAL) : sum;
        cnt[gi] = count;
        gi++;
        start = end;
    }

    return Rcpp::List::create(
        Rcpp::Named("value") = val,
        Rcpp::Named("count") = cnt,
        Rcpp::Named("key")   = key
    );
}

// ---------------------------------------------------------------------------
// Fused encode + sort + aggregate + decode in a single C++ call.
// Replaces the R-side encodeGroupKeys -> order -> grouped_agg_cpp -> decodeGroupKeys
// pipeline with zero intermediate R allocations.
//
// x:          numeric vector of values to aggregate
// group_cols: named list of grouping column vectors (each length n).
//             Columns may be integer (INTSXP) or double (REALSXP).
//             Character/factor columns must be integer-encoded by the caller.
// na_rm:      if true, NA values in x are skipped
// do_mean:    if true, compute mean; if false, compute sum
//
// Returns a named List:
//   - one element per group column (same name and same type as input),
//     containing the unique group-key values
//   - $value: numeric vector of group means/sums
//   - $count: integer vector of group sizes
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grouped_agg_from_cols(const Rcpp::NumericVector& x,
                                 const Rcpp::List& group_cols,
                                 bool na_rm = true,
                                 bool do_mean = true) {
    int n = x.size();
    int ncols = group_cols.size();
    if (ncols == 0) Rcpp::stop("group_cols must have at least one column");

    // Classify each column as integer or double and store raw pointers.
    // is_int[c] == true  -> use icols[c] (const int*)
    // is_int[c] == false -> use dcols[c] (const double*)
    std::vector<bool> is_int(ncols);
    std::vector<const int*>    icols(ncols, nullptr);
    std::vector<const double*> dcols(ncols, nullptr);

    for (int c = 0; c < ncols; c++) {
        SEXP col = group_cols[c];
        R_xlen_t clen = Rf_xlength(col);
        if (clen != n)
            Rcpp::stop("All group columns must have the same length as x");
        if (TYPEOF(col) == INTSXP) {
            is_int[c] = true;
            icols[c] = INTEGER(col);
        } else if (TYPEOF(col) == REALSXP) {
            is_int[c] = false;
            dcols[c] = REAL(col);
        } else {
            Rcpp::stop("group_cols must contain integer or numeric vectors");
        }
    }

    // Comparison lambda: compare row a vs row b across all group columns
    auto row_less = [&](int a, int b) -> bool {
        for (int c = 0; c < ncols; c++) {
            if (is_int[c]) {
                if (icols[c][a] < icols[c][b]) return true;
                if (icols[c][a] > icols[c][b]) return false;
            } else {
                if (dcols[c][a] < dcols[c][b]) return true;
                if (dcols[c][a] > dcols[c][b]) return false;
            }
        }
        return false;
    };

    auto row_eq = [&](int a, int b) -> bool {
        for (int c = 0; c < ncols; c++) {
            if (is_int[c]) {
                if (icols[c][a] != icols[c][b]) return false;
            } else {
                if (dcols[c][a] != dcols[c][b]) return false;
            }
        }
        return true;
    };

    // Build sort index.
    // For all-integer columns use a multi-pass LSD counting sort (O(n * ncols),
    // stable, cache-friendly).  Falls back to std::sort for mixed int/double.
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;

    bool all_int_cols = true;
    for (int c = 0; c < ncols; c++) if (!is_int[c]) { all_int_cols = false; break; }

    if (all_int_cols) {
        // LSD radix sort: one stable counting-sort pass per column, rightmost first.
        // pos is a stack array for small value ranges (typical SAOM data),
        // falling back to a heap vector only when range > STACK_LIMIT.
        const int STACK_LIMIT = 8192;
        int stack_pos[STACK_LIMIT];
        std::vector<int> heap_pos;
        std::vector<int> tmp(n);
        for (int c = ncols - 1; c >= 0; c--) {
            const int* col = icols[c];
            // Find value range
            int mn = col[idx[0]], mx = mn;
            for (int i = 1; i < n; i++) {
                int v = col[idx[i]];
                if (v < mn) mn = v; else if (v > mx) mx = v;
            }
            int range = mx - mn + 1;
            // Use stack array for small ranges, heap vector for large
            int* pos;
            if (range <= STACK_LIMIT) {
                pos = stack_pos;
                std::fill(pos, pos + range, 0);
            } else {
                heap_pos.assign(range, 0);
                pos = heap_pos.data();
            }
            // Count, compute start positions, scatter (stable)
            for (int i = 0; i < n; i++) pos[col[idx[i]] - mn]++;
            int acc = 0;
            for (int k = 0; k < range; k++) { int cnt = pos[k]; pos[k] = acc; acc += cnt; }
            for (int i = 0; i < n; i++) tmp[pos[col[idx[i]] - mn]++] = idx[i];
            std::swap(idx, tmp);
        }
    } else {
        std::sort(idx.begin(), idx.end(), row_less);
    }

    // Count groups
    int nG = 0;
    {
        int pos = 0;
        while (pos < n) {
            nG++;
            int cur = pos + 1;
            while (cur < n && row_eq(idx[cur], idx[pos])) cur++;
            pos = cur;
        }
    }

    // Allocate output key columns (preserve original type)
    Rcpp::List out_keys(ncols);
    // Typed pointers for output columns
    std::vector<int*>    out_icols(ncols, nullptr);
    std::vector<double*> out_dcols(ncols, nullptr);
    for (int c = 0; c < ncols; c++) {
        if (is_int[c]) {
            Rcpp::IntegerVector v(nG);
            out_icols[c] = v.begin();
            out_keys[c] = v;
        } else {
            Rcpp::NumericVector v(nG);
            out_dcols[c] = v.begin();
            out_keys[c] = v;
        }
    }
    Rcpp::NumericVector val(nG);
    Rcpp::IntegerVector cnt(nG);

    // Fill pass
    int gi = 0, pos = 0;
    while (pos < n) {
        // Store group key
        for (int c = 0; c < ncols; c++) {
            if (is_int[c])
                out_icols[c][gi] = icols[c][idx[pos]];
            else
                out_dcols[c][gi] = dcols[c][idx[pos]];
        }

        // Find group end
        int cur = pos + 1;
        while (cur < n && row_eq(idx[cur], idx[pos])) cur++;

        // Accumulate
        double sum = 0.0;
        int count = 0;
        for (int i = pos; i < cur; i++) {
            double xi = x[idx[i]];
            if (na_rm && ISNAN(xi)) continue;
            sum += xi;
            count++;
        }
        val[gi] = do_mean ? (count > 0 ? sum / count : NA_REAL) : sum;
        cnt[gi] = count;
        gi++;
        pos = cur;
    }

    // Build return list: group columns + value + count
    Rcpp::CharacterVector col_names = group_cols.names();
    Rcpp::List result(ncols + 2);
    Rcpp::CharacterVector rnames(ncols + 2);
    for (int c = 0; c < ncols; c++) {
        result[c] = out_keys[c];
        rnames[c] = col_names[c];
    }
    result[ncols] = val;
    rnames[ncols] = "value";
    result[ncols + 1] = cnt;
    rnames[ncols + 1] = "count";
    result.names() = rnames;
    return result;
}

// LOO (leave-one-out) choice probabilities for all effects in one call.
// For each effect k, zeros theta[k] and recomputes grouped softmax using:
//   util_k = util_full - contribMat.col(k) * theta[k]   (O(n) per effect)
// instead of a full matrix multiply (O(n*nEff)) each time.
// group_id must form contiguous blocks (same assumption as softmax_arma_by_group).
// Returns an nRows x nEff matrix; column k holds the LOO probs for effect k.
// [[Rcpp::export]]
arma::mat loo_change_probs(const arma::mat& contribMat,
                            const arma::vec& theta,
                            const arma::ivec& group_id) {
    int n = contribMat.n_rows, nEff = contribMat.n_cols;
    std::vector<int> starts;
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (group_id[i] != group_id[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    arma::vec util_full = contribMat * theta;  // O(n*nEff), computed once
    arma::mat out(n, nEff);
    arma::vec util(n);
    for (int k = 0; k < nEff; k++) {
        util = util_full - contribMat.col(k) * theta[k];  // O(n)
        double* out_col = out.colptr(k);
        for (int g = 0; g < nG; g++) {
            int s = starts[g], e = starts[g+1];
            softmax_inplace(util.memptr() + s, out_col + s, e - s);
        }
    }
    return out;
}

// L1 distance: sum_j |ref[j] - loo[j,k]|, aggregated by contiguous group.
// Returns nGroups x nEff; rows with <= 1 valid ref entry are set to NA.
// [[Rcpp::export]]
arma::mat l1d_grouped(const arma::vec& ref,
                       const arma::mat& loo,
                       const arma::ivec& group_id) {
    int n = ref.n_elem, nEff = loo.n_cols;
    std::vector<int> starts;
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (group_id[i] != group_id[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    arma::mat out(nG, nEff, arma::fill::zeros);
    for (int g = 0; g < nG; g++) {
        int s = starts[g], e = starts[g+1];
        int valid = 0;
        for (int r = s; r < e; r++) if (!std::isnan(ref[r])) valid++;
        if (valid <= 1) {
            for (int k = 0; k < nEff; k++) out(g, k) = NA_REAL;
            continue;
        }
        for (int k = 0; k < nEff; k++) {
            double acc = 0.0;
            const double* lc = loo.colptr(k);
            for (int r = s; r < e; r++)
                if (!std::isnan(ref[r]) && !std::isnan(lc[r]))
                    acc += std::abs(ref[r] - lc[r]);
            out(g, k) = acc;
        }
    }
    return out;
}

// KL divergence KL(ref || loo_k), normalised by log(group_size), by group.
// Returns nGroups x nEff; NA rows when group has <= 1 valid choice.
// Terms where loo_k == 0 while ref > 0 are dropped (not Inf), matching na.rm.
// [[Rcpp::export]]
arma::mat kld_grouped(const arma::vec& ref,
                       const arma::mat& loo,
                       const arma::ivec& group_id) {
    int n = ref.n_elem, nEff = loo.n_cols;
    std::vector<int> starts;
    // do we not know size beforehand or pre-allocate with a first pass?
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (group_id[i] != group_id[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    arma::mat out(nG, nEff, arma::fill::zeros);
    for (int g = 0; g < nG; g++) {
        int s = starts[g], e = starts[g+1];
        int valid = 0;
        for (int r = s; r < e; r++) if (!std::isnan(ref[r])) valid++;
        if (valid <= 1) {
            for (int k = 0; k < nEff; k++) out(g, k) = NA_REAL;
            continue;
        }
        double logN = std::log((double)valid);
        for (int k = 0; k < nEff; k++) {
            double acc = 0.0;
            const double* lc = loo.colptr(k);
            for (int r = s; r < e; r++) {
                if (std::isnan(ref[r]) || std::isnan(lc[r])) continue;
                if (ref[r] > 0.0 && lc[r] > 0.0)
                    acc += ref[r] * (std::log(ref[r]) - std::log(lc[r]));
            }
            out(g, k) = acc / logN;
        }
    }
    return out;
}

// Counterfactual mlogit update for marginal effects.
//
// Given baseline choice probabilities `p` (from softmax) and a utility
// shift `delta_u`, compute the counterfactual probability under two modes:
//
// perturbType = 0 ("alter" / one-alternative):
//   p'_j = p_j * exp(d_j) / (1 - p_j + p_j * exp(d_j))
//   This only re-weights the focal alternative j against "no change".
//
// perturbType = 1 ("ego" / ego-wide):
//   p'_j = p_j * exp(d_j) / sum_{k in group} p_k * exp(d_k)
//   This re-normalises across all alternatives in the choice set.
//
// NA values in delta_u propagate to the output as NA.
// group_id must form contiguous blocks (same assumption as softmax_arma_by_group).
// For perturbType = 0 group_id is unused but must be supplied.
//
// [[Rcpp::export]]
arma::vec mlogit_update(const arma::vec& p,
                        const arma::vec& delta_u,
                        const arma::ivec& group_id,
                        int perturbType) {
    int n = p.n_elem;
    arma::vec out(n);

    if (perturbType == 0) {
        // one-alternative update
        for (int i = 0; i < n; i++) {
            if (std::isnan(delta_u[i])) {
                out[i] = NA_REAL;
            } else {
                double ed = std::exp(delta_u[i]);
                out[i] = p[i] * ed / (1.0 - p[i] + p[i] * ed);
            }
        }
    } else {
        // ego-wide update: re-normalise within each contiguous group
        int start = 0;
        while (start < n) {
            int g = group_id[start], end = start + 1;
            while (end < n && group_id[end] == g) end++;

            // weighted = p * exp(delta_u) for this group
            double denom = 0.0;
            for (int i = start; i < end; i++) {
                double ed = std::isnan(delta_u[i]) ? 1.0 : std::exp(delta_u[i]);
                out[i] = p[i] * ed;
                denom += out[i];
            }
            if (denom > 0.0) {
                for (int i = start; i < end; i++)
                    out[i] /= denom;
            }
            // propagate NA from delta_u
            for (int i = start; i < end; i++)
                if (std::isnan(delta_u[i])) out[i] = NA_REAL;

            start = end;
        }
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector calculate_tie_prob_cpp(Rcpp::NumericVector prob,
                                           Rcpp::NumericVector density) {
    int n = prob.size();
    Rcpp::NumericVector out = Rcpp::clone(prob);
    for (int i = 0; i < n; i++) {
        double d = density[i];
        if (Rcpp::NumericVector::is_na(d) || d == 0.0) {
            out[i] = NA_REAL;
        } else if (d == -1.0) {
            out[i] = 1.0 - prob[i];
        }
        // else d == 1: keep prob[i] as-is
    }
    return out;
}