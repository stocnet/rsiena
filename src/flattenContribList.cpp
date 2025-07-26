#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame flattenContribList(List chainsList) {
  std::vector<int> chain_col, group_col, period_col, ministep_col, effect_col, choice_col;
  std::vector<std::string> effectname_col, effecttype_col;
  std::vector<double> value_col;

  int nChains = chainsList.size();
  for (int ch = 0; ch < nChains; ++ch) {
    List groupList = chainsList[ch];       // [group]
    int nGroups = groupList.size();
    for (int g = 0; g < nGroups; ++g) {
      List periodList = groupList[g];      // [period]
      int nPeriods = periodList.size();
      for (int p = 0; p < nPeriods; ++p) {
        List ministeps = periodList[p];
        int nMinisteps = ministeps.size();
        for (int m = 0; m < nMinisteps; ++m) {
          RObject mat_obj = ministeps[m];
          if (mat_obj.isNULL() || !Rf_isMatrix(mat_obj) || TYPEOF(mat_obj) != REALSXP) continue;
          NumericMatrix mat(mat_obj);
          CharacterVector effectNames = mat.attr("effectNames");
          CharacterVector effectTypes = mat.attr("effectTypes");
          int nEff = mat.nrow(), nChoice = mat.ncol();
          for (int e = 0; e < nEff; ++e) {
            std::string effname = (effectNames.size() > 0) ? as<std::string>(effectNames[e]) : "";
            std::string efftype = (effectTypes.size() > 0) ? as<std::string>(effectTypes[e]) : "";
            for (int c = 0; c < nChoice; ++c) {
              chain_col.push_back(ch+1);
              group_col.push_back(g+1);
              period_col.push_back(p+1);
              ministep_col.push_back(m+1);
              effect_col.push_back(e+1);
              choice_col.push_back(c+1);
              effectname_col.push_back(effname);
              effecttype_col.push_back(efftype);
              value_col.push_back(mat(e,c));
            }
          }
        }
      }
    }
  }
  return DataFrame::create(
    _["chain"]      = chain_col,
    _["group"]      = group_col,
    _["period"]     = period_col,
    _["ministep"]   = ministep_col,
    _["effect"]     = effect_col,
    _["choice"]     = choice_col,
    _["effectname"] = effectname_col,
    _["effecttype"] = effecttype_col,
    _["contribution"]= value_col,
    _["stringsAsFactors"]= false
  );
}