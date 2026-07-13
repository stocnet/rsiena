#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

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

// ---------------------------------------------------------------------------
// Benchmark-only softmax grouping variants — commented out 2026-07-13.
// Only softmax_rcpp_by_group (below) is used in production (predict.R).
// The variants below (ungrouped, 3-vector, matrix, list, dataframe grouping)
// were experiments benchmarked in tests/benchmarks/benchmark_rcpp.r to pick
// the fastest grouping API. Kept here (commented) so the benchmark or a future
// production need can restore them easily. To restore: un-comment the block AND
// re-add its entry to src/init.cpp, then run Rcpp::compileAttributes().
// ---------------------------------------------------------------------------
// // [[Rcpp::export]]
// Rcpp::NumericVector softmax_rcpp(const Rcpp::NumericVector& x) {
//   int n = x.size();
//   Rcpp::NumericVector out(n);
//   softmax_inplace(REAL(x), REAL(out), n);
//   return out;
// }

// [[Rcpp::export]]
Rcpp::NumericVector softmax_rcpp_by_group(const Rcpp::NumericVector& x,
                                          const Rcpp::IntegerVector& group) {
  int n = x.size();
  Rcpp::NumericVector out(n);          // one allocation for everything
  const double* xp = REAL(x);
  const int*    gp = INTEGER(group);
  double*       op = REAL(out);
  int start = 0;
  while (start < n) {
    int g = gp[start], end = start + 1;
    while (end < n && gp[end] == g) end++;              // find group boundary
    softmax_inplace(xp + start,                         // pointer into input slice
                    op + start,                         // pointer into output slice
                    end - start);                       // group size
    start = end;
  }
  return out;
}

// // [[Rcpp::export]]
// Rcpp::NumericVector softmax_rcpp_grouped(const Rcpp::NumericVector& x,
//                                          const Rcpp::IntegerVector& g1,
//                                          const Rcpp::IntegerVector& g2,
//                                          const Rcpp::IntegerVector& g3) {
//   int n = x.size();
//   Rcpp::NumericVector out(n);
//   const double* xp = REAL(x);
//   double*       op = REAL(out);
//   const int* p1 = INTEGER(g1);
//   const int* p2 = INTEGER(g2);
//   const int* p3 = INTEGER(g3);
//   int start = 0;
//   while (start < n) {
//     int end = start + 1;
//     while (end < n && p1[end] == p1[start] &&
//                       p2[end] == p2[start] &&
//                       p3[end] == p3[start]) end++;
//     softmax_inplace(xp + start, op + start, end - start);
//     start = end;
//   }
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector softmax_rcpp_grouped_mat(const Rcpp::NumericVector& x,
//                                              const Rcpp::IntegerMatrix& G) {
//   int n = x.size(), ncols = G.ncol();
//   Rcpp::NumericVector out(n);
//   const double* xp = REAL(x);
//   double*       op = REAL(out);
//   std::vector<const int*> cols(ncols);
//   for (int c = 0; c < ncols; c++) cols[c] = &G(0, c);  // contiguous per column
//   int start = 0;
//   while (start < n) {
//     int end = start + 1;
//     while (end < n) {
//       bool same = true;
//       for (int c = 0; c < ncols && same; c++)
//         if (cols[c][end] != cols[c][start]) same = false;
//       if (!same) break;
//       end++;
//     }
//     softmax_inplace(xp + start, op + start, end - start);
//     start = end;
//   }
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector softmax_rcpp_grouped_lst(const Rcpp::NumericVector& x,
//                                              const Rcpp::List& G) {
//   int n = x.size(), ncols = G.size();
//   Rcpp::NumericVector out(n);
//   const double* xp = REAL(x);
//   double*       op = REAL(out);
//   std::vector<const int*> cols(ncols);
//   for (int c = 0; c < ncols; c++) {
//     Rcpp::IntegerVector v = G[c];   // zero-copy wrap of existing R SEXP data
//     cols[c] = v.begin();
//   }
//   int start = 0;
//   while (start < n) {
//     int end = start + 1;
//     while (end < n) {
//       bool same = true;
//       for (int c = 0; c < ncols && same; c++)
//         if (cols[c][end] != cols[c][start]) same = false;
//       if (!same) break;
//       end++;
//     }
//     softmax_inplace(xp + start, op + start, end - start);
//     start = end;
//   }
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::NumericVector softmax_rcpp_grouped_cols(const Rcpp::DataFrame& data,
//                                    const std::string& val_col,
//                                    const Rcpp::CharacterVector& group_cols) {
//   Rcpp::NumericVector x_r = data[val_col];
//   const double* x = x_r.begin();
//   int n = x_r.size(), ncols = group_cols.size();
//   Rcpp::NumericVector out(n);
//   double* op = REAL(out);
//   std::vector<const int*> cols(ncols);
//   for (int c = 0; c < ncols; c++) {
//     Rcpp::IntegerVector v = data[Rcpp::as<std::string>(group_cols[c])];
//     cols[c] = v.begin();
//   }
//   int start = 0;
//   while (start < n) {
//     int end = start + 1;
//     while (end < n) {
//       bool same = true;
//       for (int c = 0; c < ncols && same; c++)
//         if (cols[c][end] != cols[c][start]) same = false;
//       if (!same) break;
//       end++;
//     }
//     softmax_inplace(x + start, op + start, end - start);
//     start = end;
//   }
//   return out;
// }

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
// build_scatter_idx — given a pre-sorted group-key matrix G_sorted (nRows × K)
// and the permutation `ord` that produced this sorted order from the original
// data, build an integer label vector of length nRows where
//   result[ord[sorted_pos]] = group_index_of(sorted_pos)
//
// Group indices are 0-based and contiguous (groups are contiguous runs in
// G_sorted, as guaranteed by grouped_agg_cpp's input convention).
//
// Used by buildAggCache to precompute row→group mappings for scatter-aggregate.
// [[Rcpp::export]]
Rcpp::IntegerVector build_scatter_idx(const Rcpp::IntegerMatrix& G_sorted,
                                      const Rcpp::IntegerVector& ord) {
    int n = G_sorted.nrow(), K = G_sorted.ncol();
    if (ord.size() != n)
        Rcpp::stop("nrow(G_sorted) must equal length(ord)");

    Rcpp::IntegerVector result(n);
    std::vector<const int*> cols(K);
    for (int c = 0; c < K; c++) cols[c] = &G_sorted(0, c);

    int gi = 0, start = 0;
    while (start < n) {
        int end = start + 1;
        while (end < n) {
            bool same = true;
            for (int c = 0; c < K && same; c++)
                if (cols[c][end] != cols[c][start]) same = false;
            if (!same) break;
            end++;
        }
        // 0-based R integer vector write — ord is 1-based (R index)
        for (int i = start; i < end; i++)
            result[ord[i] - 1] = gi;
        gi++;
        start = end;
    }
    return result;
}

// ---------------------------------------------------------------------------
// scatter_agg_1d — sequential scan of `vals` with scatter-writes into nGroups
// accumulators using a precomputed group-label vector `row_group`.
//
// row_group[i] must be in range [0, nGroups).
// na_rm: if true, NA values in vals are skipped.
//
// Returns a named list:
//   $sum   numeric vector of length nGroups — per-group sum
//   $count integer vector of length nGroups — per-group non-NA count
//
// Performance: O(n) sequential reads of vals + row_group, O(nGroups) random
// writes into accumulators.  When nGroups is small (e.g. < 1e5), the accumulators
// fit in L1/L2 cache, so this is dramatically faster than random-read permutation.
// [[Rcpp::export]]
Rcpp::List scatter_agg_1d(const Rcpp::NumericVector& vals,
                           const Rcpp::IntegerVector& row_group,
                           int nGroups,
                           bool na_rm = true) {
    int n = vals.size();
    if (row_group.size() != n)
        Rcpp::stop("vals and row_group must have the same length");

    Rcpp::NumericVector  sums(nGroups, 0.0);
    Rcpp::IntegerVector  cnts(nGroups, 0);
    const double* vp   = REAL(vals);
    const int*    gp   = INTEGER(row_group);

    for (int i = 0; i < n; i++) {
        double v = vp[i];
        if (na_rm && ISNAN(v)) continue;
        int g = gp[i];
        sums[g] += v;
        cnts[g]++;
    }
    return Rcpp::List::create(
        Rcpp::Named("sum")   = sums,
        Rcpp::Named("count") = cnts
    );
}

// ---------------------------------------------------------------------------
// grouped_agg_matrix_cpp — like grouped_agg_cpp but aggregates M outcome
// columns in a single sort pass.
//
// X:      n × M numeric matrix; each column is one outcome to aggregate.
// G:      n × K integer group-key matrix (pre-sorted by the caller, same as
//         grouped_agg_cpp expects).
// na_rm:  if true, NAs in each X column are excluded from sums/counts.
// do_mean: if true compute within-group mean; else compute sum.
//
// Returns a named List:
//   $values  nGroups × M numeric matrix — row i, col j = agg of X[,j] in group i
//   $counts  nGroups × M integer matrix — non-NA count per (group, outcome)
//   $key     nGroups × K integer matrix — unique group keys (same as grouped_agg_cpp)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grouped_agg_matrix_cpp(const Rcpp::NumericMatrix& X,
                                   const Rcpp::IntegerMatrix& G,
                                   bool na_rm  = true,
                                   bool do_mean = false) {
    int n = X.nrow(), M = X.ncol(), K = G.ncol();
    if (G.nrow() != n)
        Rcpp::stop("nrow(X) must equal nrow(G)");

    // Pointers into each column of G for group comparison
    std::vector<const int*> gcols(K);
    for (int c = 0; c < K; c++) gcols[c] = &G(0, c);

    // Pointers into each outcome column of X
    std::vector<const double*> xcols(M);
    for (int m = 0; m < M; m++) xcols[m] = &X(0, m);

    // Count groups (single pass — assumes G is pre-sorted contiguous blocks)
    int nG = 0;
    {
        int start = 0;
        while (start < n) {
            nG++;
            int end = start + 1;
            while (end < n) {
                bool same = true;
                for (int c = 0; c < K && same; c++)
                    if (gcols[c][end] != gcols[c][start]) same = false;
                if (!same) break;
                end++;
            }
            start = end;
        }
    }

    // Allocate output
    Rcpp::NumericMatrix values(nG, M);
    Rcpp::IntegerMatrix counts(nG, M);
    Rcpp::IntegerMatrix key(nG, K);

    // Fill pass — one group at a time, all M outcomes together
    int gi = 0, start = 0;
    while (start < n) {
        for (int c = 0; c < K; c++) key(gi, c) = gcols[c][start];

        int end = start + 1;
        while (end < n) {
            bool same = true;
            for (int c = 0; c < K && same; c++)
                if (gcols[c][end] != gcols[c][start]) same = false;
            if (!same) break;
            end++;
        }

        for (int m = 0; m < M; m++) {
            double sum = 0.0;
            int count = 0;
            for (int i = start; i < end; i++) {
                double xi = xcols[m][i];
                if (na_rm && ISNAN(xi)) continue;
                sum += xi;
                count++;
            }
            values(gi, m) = do_mean ? (count > 0 ? sum / count : NA_REAL) : sum;
            counts(gi, m) = count;
        }
        gi++;
        start = end;
    }

    return Rcpp::List::create(
        Rcpp::Named("values") = values,
        Rcpp::Named("counts") = counts,
        Rcpp::Named("key")    = key
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
// group_id must form contiguous blocks (same assumption as softmax_rcpp_by_group).
// Returns an nRows x nEff matrix; column k holds the LOO probs for effect k.
// [[Rcpp::export]]
Rcpp::NumericMatrix loo_change_probs(const Rcpp::NumericMatrix& contribMat,
                                     const Rcpp::NumericVector& theta,
                                     const Rcpp::IntegerVector& group_id) {
    int n = contribMat.nrow(), nEff = contribMat.ncol();
    const int*    gid = INTEGER(group_id);
    const double* th  = REAL(theta);
    std::vector<int> starts;
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (gid[i] != gid[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    // util_full = contribMat %*% theta  (dense mat-vec, computed once).
    // Accumulated column-by-column: each column of contribMat is contiguous
    // (column-major), so this is cache-friendly and needs no BLAS.
    std::vector<double> util_full(n, 0.0);
    for (int k = 0; k < nEff; k++) {
        const double* col = &contribMat(0, k);
        double tk = th[k];
        for (int i = 0; i < n; i++) util_full[i] += col[i] * tk;
    }

    Rcpp::NumericMatrix out(n, nEff);
    std::vector<double> util(n);
    for (int k = 0; k < nEff; k++) {
        const double* col = &contribMat(0, k);
        double tk = th[k];
        for (int i = 0; i < n; i++) util[i] = util_full[i] - col[i] * tk;  // O(n)
        double* out_col = &out(0, k);
        for (int g = 0; g < nG; g++) {
            int s = starts[g], e = starts[g+1];
            softmax_inplace(util.data() + s, out_col + s, e - s);
        }
    }
    return out;
}

// L1 distance: sum_j |ref[j] - loo[j,k]|, aggregated by contiguous group.
// Returns nGroups x nEff; rows with <= 1 valid ref entry are set to NA.
// [[Rcpp::export]]
Rcpp::NumericMatrix l1d_grouped(const Rcpp::NumericVector& ref,
                                const Rcpp::NumericMatrix& loo,
                                const Rcpp::IntegerVector& group_id) {
    int n = ref.size(), nEff = loo.ncol();
    const double* rp  = REAL(ref);
    const int*    gid = INTEGER(group_id);
    std::vector<int> starts;
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (gid[i] != gid[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    Rcpp::NumericMatrix out(nG, nEff);  // zero-initialized
    for (int g = 0; g < nG; g++) {
        int s = starts[g], e = starts[g+1];
        int valid = 0;
        for (int r = s; r < e; r++) if (!std::isnan(rp[r])) valid++;
        if (valid <= 1) {
            for (int k = 0; k < nEff; k++) out(g, k) = NA_REAL;
            continue;
        }
        for (int k = 0; k < nEff; k++) {
            double acc = 0.0;
            const double* lc = &loo(0, k);
            for (int r = s; r < e; r++)
                if (!std::isnan(rp[r]) && !std::isnan(lc[r]))
                    acc += std::abs(rp[r] - lc[r]);
            out(g, k) = acc;
        }
    }
    return out;
}

// KL divergence KL(ref || loo_k), normalised by log(group_size), by group.
// Returns nGroups x nEff; NA rows when group has <= 1 valid choice.
// Terms where loo_k == 0 while ref > 0 are dropped (not Inf), matching na.rm.
// [[Rcpp::export]]
Rcpp::NumericMatrix kld_grouped(const Rcpp::NumericVector& ref,
                                const Rcpp::NumericMatrix& loo,
                                const Rcpp::IntegerVector& group_id) {
    int n = ref.size(), nEff = loo.ncol();
    const double* rp  = REAL(ref);
    const int*    gid = INTEGER(group_id);
    std::vector<int> starts;
    // do we not know size beforehand or pre-allocate with a first pass?
    starts.push_back(0);
    for (int i = 1; i < n; i++)
        if (gid[i] != gid[i-1]) starts.push_back(i);
    starts.push_back(n);
    int nG = starts.size() - 1;

    Rcpp::NumericMatrix out(nG, nEff);  // zero-initialized
    for (int g = 0; g < nG; g++) {
        int s = starts[g], e = starts[g+1];
        int valid = 0;
        for (int r = s; r < e; r++) if (!std::isnan(rp[r])) valid++;
        if (valid <= 1) {
            for (int k = 0; k < nEff; k++) out(g, k) = NA_REAL;
            continue;
        }
        double logN = std::log((double)valid);
        for (int k = 0; k < nEff; k++) {
            double acc = 0.0;
            const double* lc = &loo(0, k);
            for (int r = s; r < e; r++) {
                if (std::isnan(rp[r]) || std::isnan(lc[r])) continue;
                if (rp[r] > 0.0 && lc[r] > 0.0)
                    acc += rp[r] * (std::log(rp[r]) - std::log(lc[r]));
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
// group_id must form contiguous blocks (same assumption as softmax_rcpp_by_group).
// For perturbType = 0 group_id is unused but must be supplied.
//
// [[Rcpp::export]]
Rcpp::NumericVector mlogit_update(const Rcpp::NumericVector& p,
                                  const Rcpp::NumericVector& delta_u,
                                  const Rcpp::IntegerVector& group_id,
                                  int perturbType) {
    int n = p.size();
    Rcpp::NumericVector out(n);
    const double* pp  = REAL(p);
    const double* du  = REAL(delta_u);
    const int*    gid = INTEGER(group_id);
    double*       op  = REAL(out);

    if (perturbType == 0) {
        // one-alternative update
        for (int i = 0; i < n; i++) {
            if (std::isnan(du[i])) {
                op[i] = NA_REAL;
            } else {
                double ed = std::exp(du[i]);
                op[i] = pp[i] * ed / (1.0 - pp[i] + pp[i] * ed);
            }
        }
    } else {
        // ego-wide update: re-normalise within each contiguous group
        int start = 0;
        while (start < n) {
            int g = gid[start], end = start + 1;
            while (end < n && gid[end] == g) end++;

            // weighted = p * exp(delta_u) for this group
            double denom = 0.0;
            for (int i = start; i < end; i++) {
                double ed = std::isnan(du[i]) ? 1.0 : std::exp(du[i]);
                op[i] = pp[i] * ed;
                denom += op[i];
            }
            if (denom > 0.0) {
                for (int i = start; i < end; i++)
                    op[i] /= denom;
            }
            // propagate NA from delta_u
            for (int i = start; i < end; i++)
                if (std::isnan(du[i])) op[i] = NA_REAL;

            start = end;
        }
    }
    return out;
}

// In-place contribMat → changeStats conversion for eval-only models.
// Sign-flips non-density columns on dissolution rows (density == -1) and
// optionally sets column names — all WITHOUT triggering R's copy-on-write,
// so the contribMat SEXP can be reused directly as csMat.
//
// mat:             REAL matrix (N × p), modified in-place.
// densityCol:      1-based column index of the density effect.
// newColNames:     character vector of length p (new column names, without
//                  "_eval" suffix).  Pass R_NilValue to skip renaming.
// Returns:         integer vector of length N with density values (±1/0).
// [[Rcpp::export]]
Rcpp::IntegerVector contribToCS_eval_inplace(Rcpp::NumericMatrix mat,
                                              int densityCol,
                                              Rcpp::Nullable<Rcpp::CharacterVector> newColNames = R_NilValue) {
    int n = mat.nrow(), p = mat.ncol();
    int dc = densityCol - 1;  // 0-based

    // 1. Extract density vector (integer ±1/0).
    Rcpp::IntegerVector density(n);
    for (int i = 0; i < n; i++) {
        double v = mat(i, dc);
        density[i] = ISNA(v) ? 0 : static_cast<int>(v);
    }

    // 2. Sign-flip non-density columns on dissolution rows.
    for (int j = 0; j < p; j++) {
        if (j == dc) continue;
        for (int i = 0; i < n; i++) {
            if (density[i] == -1) {
                mat(i, j) = -mat(i, j);
            }
        }
    }

    // 3. Set column names via C API (no R-level COW).
    if (newColNames.isNotNull()) {
        Rcpp::CharacterVector cn(newColNames);
        if (cn.size() == p) {
            SEXP m = mat;
            SEXP dn = Rf_getAttrib(m, R_DimNamesSymbol);
            if (Rf_isNull(dn) || TYPEOF(dn) != VECSXP) {
                dn = PROTECT(Rf_allocVector(VECSXP, 2));
                SET_VECTOR_ELT(dn, 0, R_NilValue);
                SET_VECTOR_ELT(dn, 1, cn);
                Rf_setAttrib(m, R_DimNamesSymbol, dn);
                UNPROTECT(1);
            } else {
                SET_VECTOR_ELT(dn, 1, cn);
            }
        }
    }

    return density;
}

// Softmax Jacobian matrix w.r.t. parameters.
//
// Given per-row choice probabilities p_j and the raw change-contribution
// matrix τ (contribMat) satisfying u_j = τ[j,·] · θ exactly, computes:
//
//   J[j, k] = ∂p_j / ∂θ_k = p_j * (τ[j,k] − τ̄_{G(j),k})
//
// where τ̄_{G(j),k} = Σ_{j'∈G(j)} p_{j'} · τ[j',k] is the
// probability-weighted group mean for parameter k.
//
// Arguments:
//   changeProb — [n] choice probabilities at θ̂ (output of softmax_rcpp_by_group)
//   contribMat — [n × K] raw change-contribution matrix (cc$contribMat)
//   group_id   — [n] contiguous group identifiers (same as used by softmax)
//
// Returns: [n × K] matrix of ∂p_j/∂θ_k values.
// group_id must form contiguous blocks (matching convention of softmax_rcpp_by_group).
// [[Rcpp::export]]
Rcpp::NumericMatrix softmax_jac_rcpp(const Rcpp::NumericVector& changeProb,
                                     const Rcpp::NumericMatrix& contribMat,
                                     const Rcpp::IntegerVector& group_id) {
    int n = changeProb.size();
    int K = contribMat.ncol();
    const double* cp  = REAL(changeProb);
    const int*    gid = INTEGER(group_id);
    Rcpp::NumericMatrix out(n, K);  // zero-initialized

    std::vector<double> tau_bar(K);
    int start = 0;
    while (start < n) {
        int g   = gid[start];
        int end = start + 1;
        while (end < n && gid[end] == g) end++;

        // tau_bar[k] = Σ_{j in group} p_j * contribMat[j, k]
        // (column-major iteration: each contribMat column is contiguous)
        std::fill(tau_bar.begin(), tau_bar.end(), 0.0);
        for (int k = 0; k < K; k++) {
            const double* col = &contribMat(0, k);
            double acc = 0.0;
            for (int i = start; i < end; i++) acc += cp[i] * col[i];
            tau_bar[k] = acc;
        }

        // out[j, k] = p_j * (contribMat[j, k] − tau_bar[k])
        for (int k = 0; k < K; k++) {
            const double* col  = &contribMat(0, k);
            double*       ocol = &out(0, k);
            double        tb   = tau_bar[k];
            for (int i = start; i < end; i++)
                ocol[i] = cp[i] * (col[i] - tb);
        }

        start = end;
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