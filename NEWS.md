# RSiena 1.7.0 (2026-11-18)

## Package

  * `DESCRIPTION`: website added to `URL`, and `Config/Needs/website` added.
  * `inst/CITATION`: added missing comma, and corrected Boitmanis' first name
  * Logo now also available in `man/figures/logo.png`, and added favicons
  * `NEWS.md` restructured
    - Each release date is in the version heading
    - Subheadings `Changes in RSiena` were removed
    - Section names were standardised (Effects, Functionality, 
      Bug corrections, Messages, Documentation, Coding, Tests, Package)
  * `.github/CONTRIBUTING.md` updated with more recent guidelines
  * Improved the Github CI workflows:
    - Pull request workflow (`prchecks.yml`) checks `README.md` corresponds to `README.Rmd`
    - Release workflow (`pushrelease.yml`) tags each release with the `DESCRIPTION`'s version, 
      instead of a patch increment of the newest tag
    - Release workflow attaches the manual to each GitHub release here:
      <https://github.com/stocnet/rsiena/releases/latest/download/RSiena_Manual.pdf>
    - Release workflow builds the website and deploys it to `gh-pages` branch
  * `README.md` is now generated from `README.Rmd`
    - It includes a short introduction to the SAOM, a table of contents, 
      an example that is run when the README is knitted (closes #66), 
      a section on documentation and help with introductory reading, 
      one installation section, a citation section that recommends `citation("RSiena")`, 
      and a section about the project and its funding.
  * New pkgdown website at <https://stocnet.github.io/rsiena/>
    - Configuration is in `pkgdown/_pkgdown.yml`
    - Built using colours taken from the RSiena logo
    - Function overview on the website groups help pages by task:
      making data, specifying, estimating, testing, interpreting,
      meta-analysis, and data sets
    - Navigation bar links to the scripts, the manual,
      related stocnet packages, and pages of the RSiena wiki
      (literature, applications, data sets, workshops, research projects),
      which can be updated without a new release of the package

## Documentation

  * `NEWS.md` restructured: each release date is in the version heading,
    the headings `Changes in RSiena` were removed, and the section names
    were standardised (Effects, Functionality, Bug corrections, Messages,
    Documentation, Coding, Tests, Package).

# RSiena 1.6.14 (2026-10-22)

## Functionality

  * Multi-group option (`make_group_rsiena`) now available for models
    with continuous behavior variables, with parameter starting values
    based on one fit to all groups.
  * The option `splitDepvars` in `set_algorithm_saom`, introduced in version 1.6.3,
    was rolled back and does not exist any more.

## Bug corrections

  * Analytic Phase 1 derivatives of the feedback statistic for continuous
    behavior variables in `updateDerivativeSde` corrected (it used the
    behavior at the end instead of at the start of the period).
  * There was an error when using multiple processors and maximum likelihood
    estimation. This was corrected by using `clusterApply` instead of
    `parRapply`  in `maxlike.c`.

## Tests

  * New test in `tests/newparallel.R` for multi-group data with a
    continuous behavior variable.
  * Simulation of continuous behavior data in `test16` of
    `tests/newparallel.R` corrected.

# RSiena 1.6.13 (2026-10-02)

## Effects

  * Effect `WXX` for two-mode networks corrected.
  * New effect `WXZX` (contributed by Harald Waxenecker).
  * New effect group `dyadFirstSecondBipartiteObjective`.
  * Effects  `outPop.c`, `outPopSqrt`, `outPopMore`, `outPopSqrtMore`,
    `outPopThreshold`, `inAct.c` and `inActSqrt` now have
    internal effect parameter 1.
  * Using `outPopMore`, `outPopSqrtMore`, `outPopThreshold`, or `inAct.c`
    with an internal effect parameter less than 1 leads to an error.

## Functionality

  * New method `print.coDyadCovar`.
  * New method `print.varDyadCovar`.

## Bug corrections

  * Give a constant dyadic covariate, if constructed from a sparse matrix,
    also the `vardims` attribute.

# RSiena 1.6.12 (2026-09-15)

## Effects

  * Effect `toU` corrected.
  * New effect group `dyadFirstBipartiteObjective`.
  * Effects `toU` and `WXX` now also available for two-mode networks.
    Effect group `dyadANetNetObjective` is empty now.
    (Retained for possible later use.)
  * Change to `make_specification` to allow empty effect groups.

## Functionality

  * Changes in `interpret_size` making the result a list of objects,
    starting with one object for each dependent variable;
    and dropping `interpret_size.sienaEffects`, i.e. the possibility
    to run `interpret_size` without a `sienaFit` object.
  * `toggleProbabilities` are returned by `interpret_size` only if
    `getChangeStatistics` is `TRUE`.
  * Changed the class name of the object created by `interpret_size`
    from `sienaRI` to `siena_sizes`.
  * New method `print.siena_sizes` giving the significances/p-values
    of the separate estimated parameters, the semi-standardized parameters,
    and the degree of certainty.
  * New method `print.summary.siena_sizes` giving the standard
    deviations of the change statistics and the expected relative importance
    of the effects.
  * Better error messages for `make_data_rsiena` in the case of
    incorrect node sets for bipartite networks.
  * Component `initial_theta` added to `sienaFit` objects.

## Bug corrections

  * Effect `toU` corrected.
  * Adherence to positivity constraints during parameter updates in
    Phase 2 corrected for the `doubleAveraging=TRUE` option.
  * Correction of use of `splitDepVars` in `initializeFRAN`.
  * Changes in `CalculateDerivative` in `phase1.r`
    (in a part used only by `multi_siena` to handle estimation
    difficulties).
  * In `transformScript`, deleted the last `s` in `as_covariate_rsienas`.
  * Error message in `as_covariate_rsiena`
    for square dyadic covariates corrected.

## Documentation

  * Explanation of effect groups `dyadBipartiteObjective`,
    `dyadFirstBipartiteObjective`, and `dyadSecondBipartiteObjective`,
    involving a dyadic covariate and a dependent two-mode network, added
    to manual.

# RSiena 1.6.11 (2026-07-16)

## Effects


## Functionality

  * New auxiliary functions `AC_Moran` and `AC_Geary`.
  * In `transformScript`, the replacement for `sienaCompositionChangeFromFile`
    was put before the replacement for `sienaCompositionChange`.
  * The component `theta` of results of function `siena`, the result of
    functions `coef`, and `estimate_onestep` have got a class `coef_siena`,
    which has its own `print` method (to allow vertical printing).
  * Function `coef` now has default `shortenNames=FALSE`.
  * Alias `coefficients` for `coef`.
  * Remark to console about output file for `siena` only when
    `silent=FALSE`.

# RSiena 1.6.10 (2026-06-14)

## Functionality

  * Undo sorting of returned data frames for `returnChains=TRUE`
    (`siena07models.cpp`).

## Bug corrections

  * Running `siena` with `returnChains = TRUE, returnDataFrame = TRUE`
    corrected.

# RSiena 1.6.9 (2026-06-10)

## Functionality

  * `print.sienaGroup` now also specifies the node sets for covariates
    and restricts printed line lengths to 64.
  * Use of `r` in `estimate_onestep` effectuated.
  * `checkEffects.R` corrected so that it will run without errors
    and can be used for batch checking (will take a long time).

## Bug corrections

  * `make_specification` corrected for a continuous dependent behavior variable
    and for two dependent networks and varying actor covariates (`effects.r`).
  * Use of `prevAns` with one-step estimator corrected (`initializeFRAN.r`).

# RSiena 1.6.8 (2026-06-06)

## Effects

  * Selection of effects for combinations of two dependent networks
    and an actor covariate improved, which leads to including
    more effects for combinations where one of the dependent networks
    is bipartite (`effects.r`, `effectsDocumentation.r`).

## Functionality

  * New function `estimate_onestep`.
  * Option `onestep` for `update_theta.sienaEffects` changed.
  * Hidden option to use the one-step estimator with `prevAns` in `siena`
    (`initializeFRAN.r`).
  * `print.sienadata` now specifies the node sets for covariates
    and restricts printed line lengths to 64.
  * More informative error message in `make_data_rsiena` if the number
    of elements of the second node set for bipartite networks is not correct.
  * In `transformScript`, the replacement for `descriptives.sienaGOF`
    was put before the replacement for `sienaGOF`.
  * Error message in `effectsDocumentation` if the set of `effectGroups`
    included is incorrect (meaningful only for checking the package).

## Coding

  * Reorganization and clarification of `effectGroups` for combinations
    of two dependent networks and an actor covariate:
    these now are `covarNetNetObjective`, `covarAWNetNetObjective`,
    `covarADNetNetObjective`, `covarBNetNetObjective`, `covarBXNetNetObjective`
    (`effects.r`, `effectsDocumentation.r`).

# RSiena 1.6.7 (2026-05-15)

## Functionality

  * Option `onestep` for `test_parameter(..., method="score", ...)`.
  * Option `onestep` for `update_theta`.
  * Name `sienaNet` added to `transformScript` as one of the names
    to be replaced.
  * Small comment changes in `initializeFRAN.r` and `maxlikec.r`.
  * `make_effects` also works for `sienadata` objects older than 1.4.10
    (due to small change in `effects.r`).

# RSiena 1.6.6 (2026-04-20)

## Coding

  * Stack handling of `getChangeContributionsList' improved in
  `siena07utilities.cpp', for compliance with CRAN.

# RSiena 1.6.5 (2026-04-08)

## Coding

  * In help files for `interpret_influence', `interpret_selection',
  `make_specification', `set_effect', `test_gof_auxiliary', and `test_time',
  keyword handling corrected.
  * In help file for `interpret_size_dynamics', part of example code disabled
  for testing due to time exceeding 5 seconds.

# RSiena 1.6.4 (2026-04-07)

## Bug corrections

  * For `test_gof`,
    a bug in the handling of auxiliary functions was fixed.

# RSiena 1.6.3 (2026-03-23)

## Effects

  * New effects `outdegMixedPop` and `indegMixedPop`.
  * `endowment`=TRUE for `XWX1`, `XWX2`, and `cl.XWX1`.

## Functionality

  * Auxiliary function `egoAlterCovarComb` for `test_gof`.

  * Option `splitDepvars` for `set_algorithm_saom`.
  * The quasi-Newton Raphson update step at the end of Phase 1 used
    matrix `dinvv`, which is the partially diagonalized version of `dinv`.
    This now uses `dinvv`, unless `splitDepvars`=-1.
  * Column names `inter1` and `inter2` of the result of `effectsDocumentation`
    changed to `covar1` and `covar2`.
  * Better error message for `set_interaction` in case `shortNames`
    is not a list.
  * In `sienaTest` objects, if the test is one-dimensional, the standard error
    of the linear combination is given as an additional component `sterror`
    and included in `print.sienaTest`.
  * Mean structural dissimilarity values for balance calculations
    dropped from `write_report` (almost never used, therefore confusing).
  * Several checks for objects being of class 'siena' extended to
    alternatively being of class 'sienadata'.
  * `newparallel.R' extended by a test using `test_gof`.
  * `importFrom` used for package `network`.
  * " " and "." changed to "_" in `fromObjectToText` (names for `coef`
    and `vcov`) and `fromObjectToText` (in `siena_table`).

## Bug corrections

  * For `transformScript`, processing of `Multipar.RSiena` added.
  * For `test_gof`,
    all arguments are transferred to the call of `sienaGOF`.
  * For `interpret_selection`,
    all arguments are transferred to the call of `selectionTable`.
  * For `interpret_influence`,
    all arguments are transferred to the call of `influenceTable`.

## Documentation

  * Manual updated to new function names and use of `autograph` for
    plotting influence and selection tables.
  * In help file for `set_interaction`, treatment of `shortNames`
    corrected.
  * In help file for `set_effect`, treatment of `covar1` and `covar2`
    corrected.

## Coding

  * `Rf_PrintValue` in `Chain.cpp` dropped
    (was leftover from earlier debugging).
  * `Rf_findVar` replaced by `R_getVar` in `siena07models.cpp`
    in view of requirements of `R 4.5.3`; requires at least `R 4.5.0`.
  * Depends on R changed to `(>= 4.5.0)`.

# RSiena 1.6.2 (2026-02-26)

## Functionality

  * Function `transformScript`.

  * New method `write_report` for `sienadata` objects,
    replacing `print01Report` (which still is retained).

  * Attribute `EgoAlter` added to `sienaGOF` results,
    which is `TRUE` currently only for `egoAlterCombi`
    (the purpose is to use it for plotting).

## Coding

  * `returnDeps` mentioned explicitly as an argument for `siena`
    (and not just as an element of ...).
  * Some occurrences of deparse(substitute(.)) accompanied by
    a check and subsequent potential stop in case of
    objects resulting from piping.
  * `sienaAlgorithmCreate` modified for different operation
    in case it is called from `transformScript`.
  * 'parallel.R' in the `tests` directory is replaced by 'newparallel.R'
    which uses the new names and no `parallelTesting`.

# RSiena 1.6.1 (2026-01-20)

## Effects

  * Effects `sameXOutAct`, `diffXOutAct`, `crossXOutAct` also implemented
    for two-mode networks.
  * Parameter values 0 and -1 added for `outPop` and `inAct`
    (earlier these were wrongly mentioned here as `outAct` and `inPop`).

## Coding

  * <math.h> replaced by <cmath> in C++ code.

# RSiena 1.6.0 (2026-01-04)

## Functionality

  * Overhaul of naming system, which required many new function and methods.

  * Warning if `egoAlterCombi` is used with more than 99 ego-alter values.

## Bug corrections

  * `checkImpossibleChanges` restricted to networks.
  * `checkZeroChanges` adapted to also cover behavioral variables.

## Documentation

  * Chapter 3 was inserted in the manual to explain the new function names.
    Appendix A, "List of Functions in Order of Execution", was dropped
    and the author of this appendix, Paulina Preciado Lopez, was
    accordingly dropped from the list of authors of the manual.
    Greetings, Paulina!}
  * Example of `dyadicCov` for ego-alter combinations dropped,
    as this is already covered by `egoAlterCombi`.

# RSiena 1.5.7 (2025-11-03)

## Effects

  * Completely refactored diffusion rate effects as a child of newly
    constructed abstract class behaviorRateEffect that is used in
    StatisticCalculator and DependentVariable for calculation of rates, scores
    and statistics. New diffusion rate effects can now be added relatively
    easily with a new cpp file. Closes #97 .
  * NetworkDependentBehaviorEffects can now also use twopath tables.
    Also added a two-path-count helper function.
  * New distance 2 exposure diffusion rate effects `totInExposureDist2`,
   `totAInExposureDist2`, `anyInExposureDist2`, behavior effects `totPopAlt`,
   `totGroup`, `indegAvGroup`, `indegtotGroup`, `totGwdspFFAlt`, `totGwdspFBAlt_nc`
   and new behavior effects without centering `quad_nc`,
   `avInAltDist2_nc`, `totInAltDist2_nc`, `avTinaltDist2_nc`, `totAInAltDist2_nc`,
   `totGwdspFFAlt_nc`, `totGwdspFBAlt_nc`.
 * New effects `avSameXAlt`, `totSameXAlt`, `avSameXInAlt`, `totSameXInalt`,
   `avXRecAlt`, totXRecAlt`, `avSameXRecAlt`, `totSameXRecAlt`,
   `avAltSameX`, `totAltSameX`, `avRecAltSameX`, `totRecAltSameX`.
  * New effects `avInAltAltX`, `totInAltAltX`, `avRecAltAltX`, `totRecAltAltX`,
    `avAltSameX`, `totAltSameX`, `avInAltSameX`, `totInAltSameX`,
    `avRecAltSameX`, `totRecAltSameX`.

## Documentation

  * A new subsection was added to `Siena_algorithms` to explain the new
    construction of diffusion rate effects.

## Coding

  * Slightly refactored contribution, target and actortarget getter to use
  a unified API for setting up C++ objects.
  * Corrected how static contributions are extracted to actually calculate
    correct contributions for the first potential ministep in each period,
    using the values at the start of the period to calculate allowed changes and
    effect change statistics. This means that there is no calculation of the
    first ministep after the last observation wave anymore. This also means that
    sienaRI does not calculate a RI measure at the last observation moment.
    Added clarification in and cleaned up the documentation and scripts
    for `siena07`, `sienaRI`, and `sienaRIDynamics`.
  * Eliminated additional test file for chains and integrated into `parallel.R`
    as far as it seemed necessary.
  * Cleaned up some legacy code related to rate effects.
  * Reprogrammed effects `avXAlt`, `totXAlt`, `avXInAlt`, `totXInalt`
    for greater efficiency (class `AltersCovariateAverageEffect`).
    Removed class `InAltersCovariateAverageEffect`.
  * Dropped `EpochSimulation::updateParameters`,
    `DependentVariable::updateBasicRate`,
    `DependentVariable::updateEffectParameters`. These were not used at all.
  * Dropped algorithm as an argument for function `unpackData` and
    `getChangeContributions`, and as a consequence also from
    `getTheActorStatistics` and argument `x` from `actorTargets`.

# RSiena 1.5.6 (2025-10-28)

## Functionality

  * Auxiliary function `egoAlterCombi` for GoF now omits values
    where the behavioral variable is `NA`.
  * Various small additions to output of `siena07` for `verbose=TRUE`
    (`phase1.r` and `phase2.r`).

## Coding

  * The virtual function `preprocessEgo` now is defined properly for
    the function classes `NetworkAlterFunction`,
   `CovariateNetworkAlterFunction`, `DoubleCovariateFunction`,
   `SameCovariateInTiesFunction`.
   This was done to correct a bug that appeared for effect `sameXInPop`
   for a two-mode network; perhaps there was a similar bug for `diffXInPop`,
   `sameXOutAct`, `diffXOutAct`, and `sameXVInPop`.

# RSiena 1.5.5

## Functionality

  * Auxiliary function `egoAlterCombi` for GoF now omits values
    where the behavioral variable is `NA`.
  * Various small additions to output of `siena07` for `verbose=TRUE`
    (`phase1.r` and `phase2.r`).

## Coding

  * The virtual function `preprocessEgo` now is defined properly for
    the function classes `NetworkAlterFunction`,
   `CovariateNetworkAlterFunction`, `DoubleCovariateFunction`,
   `SameCovariateInTiesFunction`.
   This was done to correct a bug that appeared for effect `sameXInPop`
   for a two-mode network; perhaps there was a similar bug for `diffXInPop`,
   `sameXOutAct`, `diffXOutAct`, and `sameXVInPop`.

# RSiena 1.5.3 (2025-10-17)

## Effects

  * New effects `quad_cc`, `avAlt_cc`, `totAlt_cc`.

# RSiena 1.5.2 (2025-09-06)

## Effects

  * New effects `outActMore_ego`, `outActSqrtMore_ego`, `outMore_ego`,
    `outPopMore`, `outPopSqrtMore`, `outPopThreshold`.

## Functionality

  * `selectionTable` got an attribute `quad` indicating whether the plot
    is a quadratic function.
  * Users can now extract now extract the changeContributions when running
    `siena07` by setting an argument `returnChangeContributions=TRUE`.
    If used together with `nsub=0` and `prevAns` or modified initial
    values in the effects object, especially useful for post-estimation,
    e.g. in `sienaRIDynamics`.
  * `sienaRIDynamics` uses siena07 directly now and is reinstated.

## Coding

  * Use `(any(!gmm))` in  `initializeFRAN` to allow the use of a `prevAns`
    object with a different estimation method for a multigroup estimation.

# RSiena 1.5.1 (2025-07-12)

## Effects

  * New effects `divOut_ego`, `divIn_ego`.

## Functionality

  * `sienaRI` reinstated.
  * Option `prML=2` for maximum likelihood estimation
    using the `move` proposal step reinstated (`sienaAlgorithmCreate`).

## Coding

  * Imported functions from packages `Matrix`, `lattice`, `parallel`, `MASS`,
    and `xtable` mentioned specifically in the `Namespace`
    instead of importing these entire packages.

# RSiena 1.5.0 (2025-07-05)

*New CRAN version.*

## Package

  * Christian Steglich now is maintainer.

# RSiena 1.4.25 (2025-07-05)

## Coding

  * Undid other changes in version 1.4.23 to
    `getTargetsChangeContributions` in `siena07setup.cpp`.

# RSiena 1.4.24 (2025-07-04)

## Effects

  * Error for creating effects in effects group `dyadSecondBipartiteObjective`
    was corrected in file `effects.r`,
    so they now are included also for changing dyadic covariates.
    This affected effects `XWX`, `XWX1`, and `XWX2`.

## Coding

  * Function `Chains:printConsecutiveCancelingPairs` deleted from `model\ml`,
    because it was not used and led to a protection error.
  * Undid changes in version 1.4.23 to `getTargetActorStatistics` and
    `getTargetsChangeContributions` in `siena07setup.cpp`.

# RSiena 1.4.23 (2025-05-03)

## Effects

  * New effects `fromAny`, `sameInXCycle4`.
  * Effect shortName `cycle4ND` replaced by `cycle4`.
  * Better treatment of missing covariate values in effect `sameXCycle4`.
  * Internal effect parameter values 3 and 4 for `sameXInPop`,`diffXInPop`,
    `sameXInPopIntn`,  `sameXInActIntn`, `homXOutAct2`.
  * New C++ class `CatCovariateDependentNetworkEffect` (for `homXOutAct2`).
  * Added `(#)` to the `effectName` of `from`, `sameXInPop`, `diffXInPop`,
	`sameXVInPop`, `sameXVInPop2`.

## Functionality

  * In `print.sienaEffects(..., includeShortNames=TRUE)`, the `effectNumber`
    is also printed.

## Bug corrections

  * Repaired (hopefully) memory leak in `getTargetActorStatistics` and
    `getTargetsChangeContributions` in `siena07setup.cpp`.
  * Repaired (hopefully) memory leak in `move` in `MLSimulation.cpp`.
  * Repaired memory leak  in `DoubleCovariateCatFunction`.

## Documentation

  * Definition of effects `outRateLog`, `inRateLog`, and `recipRateLog`
    corrected in the manual.

## Coding

  * Class `Covariate` has new variable `covariateN`, which then is transferred
    as `covarN` to `CovariateDependentNetworkEffect` and
    `DoubleCovariateFunction`.

# RSiena 1.4.22 (2025-02-08)

## Effects

  * Internal parameter 0 (for log(x)) added for `outActIntn`.
    (It would be trivial to implement this also for the other
    mixed degree effects, but currently there seems no need.)

## Functionality

  * New auxiliary GOF function `egoAlterCombi`.
  * Parameter `showAll` added to `plot.sienaGOF`.

  * New table `IntLogTable` and new generic function`IntLogFunction`.

## Bug corrections

  * In `sienaGOF`,  if the auxiliaryFunction does not always
    give vectors of the same length, the error message gives properly
    the name of the auxiliaryFunction.

# RSiena 1.4.21 (2024-12-18)

## Effects

  * New effect `altHigherEgoX`.
  * Effect `higher` also implemented for symmetric networks.
  * Used # in name of threshold effects.

## Functionality

  * New functions `selectionTable` and `influenceTable`,
    with `print` methods.

## Bug corrections

  * In `siena.table`,  fixed parameter values are not reported as NA,
    but as their fixed values.

## Messages

  * If `includeInteraction` is called with argument `parameter`,
    an error message appears that this keyword should not be given.
  * Extension of message in `sienaTimeTest` given in the case of
    collinearities.

## Documentation

  * Improved explanation of the use of `interaction1` and `interaction2`
    in the help page for `includeInteraction`.

# RSiena 1.4.20 (2024-11-10)

## Effects

  * New effects `varAlt` and `avSimVarAlt`.

## Bug corrections

  * Interaction effects for continuous behavior now turned off. The R side of
    `includeInteraction` for continuous behavior had been implemented but
    not the C++ side, so now an error message is given.

## Messages

  * If `includeInteraction' is called for a continuous behavior effect, an
    error message is given that interaction effects are not yet implemented
    for continuous behavior.

## Coding

  * For the new effects, `lvariance`, the behavior variance over all but the
    last wave, was included as member data of `BehaviorLongitudinalData`
    together with the corresponding member functions.

# RSiena 1.4.19 (2024-09-03)

## Effects

  * New effects `outThreshold` and `outThreshold2`.
  * `outActIntn` is an ego effect (`allEffects.csv`).

## Functionality

  * Network option 11 = `NETCONTEMP` for use of contemporaneous statistics
    for estimating all evaluation effects of the network variable.

## Bug corrections

  * Error corrected that occurred for estimation by the GMoM if
    some effects are fixed.
  * List of effect names when `thetaBound` is exceeded corrected (`phase2.r`).

## Messages

  * If `setEffect` is called for an effect with `type=gmm`, an error message
    is given that the function called should be `includeGMoMStatistics`.

# RSiena 1.4.18 (2024-07-31)

## Effects

  * New gmm effects for co-evolution of multiple networks:
    `crprod_gmm`, `to_gmm`, `from_gmm`.

## Bug corrections

  * Bug corrected in `sienaGOF` which occurred for models with tested effects
    if `iterations` is less than `sienaFitObject$n3`.

## Coding

  * The new gmm effects required the extension of the function classes
    `AlterFunction`, `NetworkAlterFunction`, and `MixedNetworkAlterFunction`
    to allow making estimation statistics depend totally on the simulated state.
  * Function `fixUpEffectNames(effects, defaultEffects)`
    moved from `initializeFRAN.r` to `sienaEffects.r`.

# RSiena 1.4.13 (2024-06-06)

## Effects

  * New effect `homXOutAct2`.

## Functionality

  * The handling of internal effect parameters in the columns `effectName`
    and `functionName` of `sienaEffects` objects is improved.
    This was achieved by changes in functions `setEffect`
    and `includeInteraction`.
    Effectnames as reproduced by `sienaFit` objects should now contain the
    correct values of internal effect parameters.
  * New parameter `thetaBound` for `siena07`, which has the effect of stopping
    the estimation process if some parameters become too large (which would
    signal divergence).
  * Display of deviations in `siena07gui` modified so that numbers larger
    than or equal to 1e5 in absolute value are displayed in exponential format
    (and use only one line in the gui) (function `FormatString` in `siena07.r`).

## Bug corrections

  * Function `setEffect` corrected (it did not give the proper internal
    parameters in the `effectName`).
  * Function `updateSpecification` corrected (it did not work properly
    for including interactions).

# RSiena 1.4.12 (2024-04-27)

## Coding

  * All objects created by functions, if not print or summary,
   now have an attribute "version", which is the package version.

## Package

  * New RSiena logo inst/rsienalogo-2.png

# RSiena 1.4.11 (2024-04-23)

## Effects

  * Effect name `outOutDist2AvIntn` changed to `avAlt.2M.tot`.
  * New effects `avAlt.2M.tie`, `avAlt.2M.tot`, `avAltU.2M.tie`,
    `dist2OutInActIntn`, `nDist2ActIntn`, `sharedToU`.
  * Changed sqrt treatment in `outOutDist2ActIntn` and
    `outOutDist2AvIntn` / `avAlt.2M.tot`.

## Functionality

  * Method `print.sienaEffects` has an extra parameter `includeShortNames`
    to do what the name of this parameter suggests.
  * Improved error message for function `updateSpecification`.

# RSiena 1.4.10 (2024-03-25)

## Effects

  * Effects `sameXV` and `sameXVInPop` added for symmetric networks,
    and restricted to integer-valued variables in the range from 0 to 20.

## Functionality

  * Parameter `silent` (new in version 1.4.8)
    in `sienaAlgorithmCreate` activated.

## Coding

  * Actor covariates in `sienaData` have a new attribute `lowIntegers`
    used for in/excluding effects `sameXVInPop` and `sameXV` in
    `getEffects`.

# RSiena 1.4.9 (2024-03-21)

## Effects

  * New effects `crossXOutAct`, `outOutDist2ActIntn`,
    `outOutDist2AvIntn`, `inPopOutW`.
  * New effect group `doubleCovarNetObjective`.
  * New effects `sameXV` and `sameXVInPop` for bipartite networks.
  * `sameXCycle4` added for one-mode and symmetric networks.
  * `sharedTo` gets default internal effect parameter `p=3`.

## Functionality

  * Function `updateSpecification` now also updates interaction effects
    and `initialValues`.

## Coding

  * Internal functions `numberIntn`, `numberBehIntn`, `checkVersion`
    defined in file `initializeFRAN.r`.
  * The `Covariate` class and its descendants (all actor covariates)
    now have functions `min` and `max`.

# RSiena 1.4.8 (2024-02-29)

## Functionality

  * New parameter `silent` in `sienaAlgorithmCreate`.

## Bug corrections

  * Correction of memory leak in `siena07setup.ccp` for ML estimation.

# RSiena 1.4.7 (2024-02-21)

*New CRAN version.*

## Bug corrections

  * Correction (again) of error message in `siena07utilities::Rterminate`.

# RSiena 1.4.6 (2024-02-19)

## Functionality

  * New parameter `targets` in `siena07`, used to supersede the targets
    calculated from the data (not for use in estimation for regular data sets,
    see the help file for `siena07`).
  * `effectsDocumentation` reports to the console
    the name of the file that was written.
  * `sienaDataCreate` stops with an error message if there is a bipartite
    network before a one-mode network.
  * If the data contains a continuous dependent behavioral variable and
    the algorithm specifies conditional estimation, the estimation stops
    with a clear error message.

  * `sienaRI` temporarily disabled because of a memory leak.
  * The option `prML=2` temporarily disabled because of a memory leak in the
    `move` proposal distribution (`sienaAlgorithmCreate`).
  * `doMoreUpdates` moved from `maxlike.r` to `maxlikecalc.r`.
    The rest of `maxlike.r` as wel as `maxlikefn.Rd` deleted.
    These were not used anywhere.

## Bug corrections

  * R_NO_REMAP included, and "\Rf_" prepended to all function names used
    from `Rinternals.h` and `R_ext/Error.h` (in `siena07models.cpp`,
    `siena07utilities.cpp`, `siena07setup.cpp` and `siena07internals.cpp`)
    and various other places
    (in accordance with "Writing R Extensions").
  * Corrected wrong length of `lprobabilityArray` in `MLSimulation.h`,
    and cleaned up a bit.
  * In `MLSimulation.cpp`, various sets of "delete" commands reordered so as
    to be in opposite order of the corresponding "new" commands.
  * Correction of error message in `BehaviorVariable::accumulateDerivatives`
    and in `siena07utilities::Rterminate`.
  * Allow `sienaDataCreate` to work with a single variable defined as
    a dependent network given as a list of sparse matrices.
  * Allow `getEffects` to construct effects of more than one dependent network
    on continuous behavior dependent variables.
  * Added some "drop=FALSE" in `initializeFRAN.r` to guard against
    dimension loss in the construction of sparse matrices.
  * Change check for constant dyadic covariates for sparse matrices (issue #88).
  * Some of the recently added effect groups were missing from
    `effectsDocumentation`. This led to an incomplete listing of the effects.
    They are now included.

## Package

  * List of changes in the code before 2022 moved from `NEWS.md` to `ONEWS_gh`.

# RSiena 1.4.5 (2024-02-14)

## Bug corrections

  * Put `#include <Rinternals.h>` as the last of the include commands
    in various `.cpp` files (Tomas Kalibera).

## Package

  * Migrated package repository to "stocnet" organisation.

# RSiena 1.4.2 (2023-12-14)

## Effects

  * New effects `outXMore`, `outMore3`.
  * `Interactiontype` of `altLThresholdX` and `altRThresholdX` is dyadic.
  * `Interactiontype` of `degAbsDiffX`, `degPosDiffX`, and `degNegDiffX`
    is "" (blank).
  * Corrected effect `outMore`.

## Functionality

  * New parameter `iterations` in `sienaGOF` to allow shorter computations.

## Coding

  * `bxeffects` initialized to 0 in `ContinuousVariable::accumulateScores`.
  * All private variable declarations in the C++ `*.h` files
    were initialized using braces.
  * In `mixedTriadCensus`, the check of the condition about the matrix
    dimensions was split into its two parts.
  * `TruncatedOutdegreeEffect2` dropped from `src/model/effects`
    (was superfluous).

# RSiena 1.4.1 (2023-11-01)

*New CRAN version.*

## Documentation

  * Changes in accordance with "Guidelines for Rd files".

# RSiena 1.4.0 (2023-10-31)

## Coding

  * Superfluous "lsimulatedDistance" in `SdeSimulation.h` deleted.
  * In `PrimarySettingEffect.cpp`, used `to_string` for converting a number
    to string in an error message.

## Package

  * `Author` field omitted in `DESCRIPTION`, because `Author@R` is sufficient.

  * Superfluous files in `inst` directory omitted.
  * New `CITATION`.

# RSiena 1.3.28 (2023-10-11)

## Effects

  * `threshold`, `threshold2`, `threshold3`, `threshold4` changed to
    work with non-centered parameters (not backward-compatible).

## Documentation

  * Descriptions of effects `altInDist2W` and `totInDist2W` added
    to the manual (the effects had been there since a long time,
    but not documented).

# RSiena 1.3.27 (2023-09-29)

## Effects

  * `avInSimDist2`, `totInSimDist2`, `sameEgoDist2`,`sameEgoInDist2`,
    `outMore2`, `divOutEgoIntn`, `divInEgoIntn`, `divOutAltIntn`,
    `divInAltIntn`.
  * `avTAltDist2` and `totTAltDist2` also implemented for behavior co-evolving
    with symmetric networks.

## Documentation

  * Some explanation is given in the manual about internal effect parameters
    for interactions created by `includeInteraction`.

## Coding

  * Corrected one line in `siena07models.cpp`, which led to slowness of
    `siena07` since version 1.3.18.

# RSiena 1.3.26 (2023-08-15)

## Bug corrections

  * Period/groupwise tests in `sienaTimeTest` corrected for the case of
    non-saturated sets of dummy variables.
  * `plot.sienaTimeTest` for "pairwise=TRUE" changed so that the warning
    is avoided.
  * `sienaGOF` corrected so that again it can handle auxiliary functions
    referring to more than one `varName` (such as in `mixedTriadCensus`).

## Coding

  * Improved Phase 1 derivative matrix computation for basic SDE parameters.
  * Added continuous behavior to returned simulated data.

# RSiena 1.3.24 (2023-08-01)

## Bug corrections

  * In `getEffects`, the effects object was given an attribute `version`,
    which was not done correctly in version 1.3.23.
    (This led to always giving a warning if any interaction effects were
    specified.)
  * Corrections of implementation of acceptance by `sienaGOF` of a list of
    `sienaFit` objects (was not correct in version 1.3.23).

## Tests

  * function `includeInteraction` used in "parallel.R" (for testing).

# RSiena 1.3.23 (2023-06-29)

## Effects

  * New effects `diffWXClosure`, `sameWWClosure`,  `diffWWClosure`,
    `diffXWClosure`, `sameXWClosure`, `unequalX`.
  * `JoutMix` made available for bipartite dependent networks.
  * For continuous behavior variables depending on a bipartite dependent
    network, the effect group `continuousBipartiteObjective` was created,
    with effects `outdeg`, `outdegSqrt`, and `isolateOut`.
  * `sameXOutAct` and `diffXOutAct` now have a parameter 2 for `sqrt`.

## Functionality

  * The model for continuous behavior variables seems to work now,
    because of the first correction mentioned above.
  * `sienaGOF` now also accepts a list of `sienaFit` objects.

## Bug corrections

  * In `initializeFRAN`, the call of `getEffects` now is dependent on
    the value of attribute "onePeriodSde".
  * The error was corrected that occurred if `useStdInits = TRUE`
    in `sienaAlgorithmCreate` and the effects object includes
    interaction effects.
  * In `sienaDataCreate`, the warning message that there is at least one
    `upOnly` period now is made for each dependent variable instead of
    only the last.
  * In `getEffects`, the effects object was given an attribute `onePeriodSde`
    and an attribute `version`.
  * In `initializeFRAN`, the comparison between `effects` and
    `defaultEffects` now is based on `shortName` instead of `effectName`
    (`effectName` was changed if there are interaction effects),
    excluding the lines in the effects object for `unspInt` and `behUnspInt`
    to allow effects objects created
    with non-default values of `nintn` and `behNintn`.

## Coding

  * Better text for stop in `initializeFRAN` when there is a mismatch
    between effects objects disabling the creation of interaction effects.
  * Warning in `initializeFRAN` if the version of the effects object
    is not current and the effects object contains interaction effects
    (then it is possible that the interacting effects are chosen incorrectly,
    even though the `effectName` of the interaction seems OK).
  * Better error message in `sienaGOF` if `groupName` or `varName` is incorrect.
  * Use default bandwidth selection in violin plot for `sienaGOF`
    (the use of "nrd" sometimes led to absent plots because of negative bw).

# RSiena 1.3.22 (2023-05-11)

## Functionality

  * For one-mode networks, new model options `DOUBLESTEP25`, `DOUBLESTEP50`,
    `DOUBLESTEP75`, `DOUBLESTEP100`.

## Bug corrections

  * The first item in "Coding" implies correction of several distance-2 network
	 effects such as `altDist2`, `totDist2` and `altInDist2`.
  * In `sienaAlgorithmCreate`, changed default `prML=2` back to `prML=1`;
    stop if Maximum Likelihood estimation is attempted for a data set
    containing more than one dependent variable
    with `prML=2` (implemented in `initializeFRAN.r`).

## Coding

  * Corrected and cleaned up virtual definitions in `AlterFunction` and its
    descendants, in particular `CovariateNetworkAlterFunction`.
  * Added `const` to virtual specification of `value` in `AlterFunction.h`
    and all of its descendants.
  * Replaced ambiguous call to `std::abs` in `AbsDiffFunction.cpp`.

# RSiena 1.3.20 (2023-04-22)

## Functionality

  * The keyword `parameter` in `includeInteraction` was dropped because it did
    not have any consequences. The help page for `includeInteraction` now
    explains how internal effect parameters for user-defined interactions
    are determined.
  * The column `dimnames` of the `Simulations` array returned by `sienaGOF`
    are set to the names of the elements of the auxiliary function.
  * Standard deviations added to output of `descriptives.sienaGOF`.
  * Improved error message in `initializeFRAN` in the case of mismatch between
    effects objects.
  * Warning in `sienaAlgorithmCreate` if `(maxlike && (!is.null(MaxDegree)))`.
    This is now also mentioned in the help page for `sienaAlgorithmCreate`.

## Bug corrections

  * `updateSpecification` (in `effectsMethods`) now also updates
    internal parameter values.
  * In `TriadCensus`, the empty network will not lead to an error
    but be reported with the correct triad census.
  * For `reciAct`, check whether internal parameter ==2 replaced by check
    whether absolute difference from 2 is less than 0.001.
  * In `phase2.r`, `z$sd` is calculated using `sqrt(pmax(..., 0))` to avoid the
    extremely rare case of a negative calculated variance.
  * In `sienaDataCreate`, handling of structurally determined values
    in `checkConstraints` corrected (thanks to issue raised by Jos Elkink).

## Documentation

  * Reference about score-type test added to `Wald.Rd`.
  * In the help page for `sienaDependent`, it is mentioned that if there are
    one-mode as well as two-mode dependent networks,
    the one-mode networks should come first.

# RSiena 1.3.19 (2023-02-07)

## Effects

  * New effect `inPop_dya`.
  * Parameter 2 for `sameXInPop` and `diffXInPop`.

## Functionality

  * `coCovar` and `varCovar` now can handle variables with only one
    non-missing value, but will stop with an error message
    if all values are missing.

## Bug corrections

  * Help page for `siena07` corrected with respect to `x$lessMem`.

## Coding

  * `siena07internals.cpp` adapted to be compatible with new clang 16 C++
    compiler (thanks to Brian Ripley).

# RSiena 1.3.18 (2023-01-29)

## Functionality

  * Additional step type `move` for MH proposal distribution
    for likelihood estimation (thanks to Charlotte Greenan).
  * Accordingly, parameters changed that are used in `sienaAlgorithmCreate`
   for probabilities of MH steps, now summarized in `prML`; with a new default.
  * List elements `accepts`, `rejects`, `aborts` for `sienaFit` objects
    produced by ML estimation improved/corrected by reorganizing them in C++.
  * List element `ac3` added to `sienaFit` object if `maxlike`.

# RSiena 1.3.17 (2023-01-06)

## Functionality

  * `sienaGOF` now accepts simulated auxiliary statistics containing missing
    values. If there are any, this will be reported with a warning
    if `giveNAWarning` is `TRUE`.
  * `sienaDataCreate` now also accepts, as "...", a list of such objects.

# RSiena 1.3.16 (2023-01-02)

## Effects

    `inPopIntnX`, `inActIntnX`, `outPopIntnX`, `outActIntnX`, `sameXInPopIntn`,
    `sameXOutPopIntn`, `sameXInActIntn`, `sameXOutActIntn` restored
    (these had got lost in some way...).

## Bug corrections


## Coding

  * All occurrences of `http` in `R` and `Rd` files changed to `https`.
  * `seq_len` used and superfluous `c()` omitted in various R files.

# RSiena 1.3.15 (2022-11-27)

## Bug corrections

  * `siena08`: correct p-value `pTsq` for overall test statistic `Tsq`
  * `print.summary.sienaMeta`, `siena07`, `print01Report`: drop RForge revision.
  * Correct "objname" to "projname" in `meta.table` (`siena08.r`).
  * Simplify `LaTeX` output of `meta.table`.
  * `seq_along` and `seq_len` used in `print01Report`.

# RSiena 1.3.14 (2022-11-04)

*New CRAN version.*

## Bug corrections

  * Update `configure` and  `configure.ac` (with help from Brian Ripley).

# RSiena 1.3.13 (2022-10-07)

## Coding

  * Replacements in EffectFactory.cpp of single | operator by ||.

# RSiena 1.3.12 (2022-10-06)

## Effects

  * Several new effects related to primary setting:
    `nonPCompress`, `primCompress`, `primary`, `primDegAct`,
    `primDegActDiff`, `primDegActDiffSqrt`, `primDegActSqrt`,
    `primDegActLog`, `primDegActInv`.
  * `gwdspFB` effect added for two-mode networks.
  * New effects `outAct_ego`, `inAct_ego`,`reciAct_ego`, `toAny`.
  * For effects `to`, `toBack`, `toRecip`, `mixedInXW`,
    internal effect parameter 3 now specifies truncation of the number of
    twosteps (change to `MixedTwoStepFunction`).

## Bug corrections

  * `universalOffset` initialized as 0; it was earlier initialized as
    the maximum real number (`NetworkLongitudinalData.cpp`).
  * `thetaStore` deleted (was trash in `phase2.r`).
  * Various comparisons for vectors with 0 changed to using `all`
    to avoid warnings (`initializeFRAN.r`).

## Documentation

  * Modified help page for `sienaRI`.
  * Small modifications of help page for `sienaGOF`.

## Coding

  * Changes to comply with new version of `Matrix` package.
  * Replacements in some C++ functions of single & and | operators by && and ||.

  * `sigmas` and `meansigmas` added to `sienaRI` object.
  * Print of standard deviations in the `sienaRI` object for `printSigma=TRUE`
    changed to using averages at the variance level.
  * If `returnThetas` in the call of `siena07`, also simulated estimation
    statistics during Phase 2 (deviations from targets) are returned.

# RSiena 1.3.11 (2022-05-30)

## Functionality

  * Further explanation of `mixedTriadCensus` in the help page for
    `sienaGOF-auxiliary`.

## Bug corrections

  * Correction in `effects.r` of error that led to warning
    for multivariate networks.
  * Correction of help page for `sienaGOF` (`groupName`).
  * Correction of `igraphNetworkExtraction` in the help page for
    `sienaGOF-auxiliary`.

# RSiena 1.3.10 (2022-04-28)

## Effects

  * New effects `avInAltW`, `avWInAlt`, `totInAltW`, `totWInAlt`
    (with help from Robert Krause).
  * Corrected implementation of `sharedTo`.

## Bug corrections

  * Bug corrected that occurred when several two-mode networks were included
    in the dependent variables, with an order restriction between them.
    (Correction of `HigherFilter` and `DisjointFilter`).

## Coding

  * Several modifications to enable traceback of errors occurring
    in checkSenderRange called in inTies.

# RSiena 1.3.9 (2022-03-18)

## Effects

  * Corrected implementation of `simAllNear` and `simAllFar`.

## Bug corrections

  * small correction of `summary.sienaGOF`.
  * small correction of `sienaTimeTest`.

# RSiena 1.3.8 (2022-03-07)

## Effects

  * Changed default internal effect parameter for `simAllNear` to 2 and for
    `simAllFar` to 4.

## Functionality

  * In `sienaTimeTest`, added `warn=FALSE` to `varCovar()` to avoid warnings.
  * Small improvements to help pages for `sienaGroupCreate` and `sienaGOF`.

## Bug corrections

  * corrected sienaRI for behavioral variables.
    This required changes in
    `StatisticCalculator::calculateBehaviorStatistics` and
    `StatisticCalculator::calculateBehaviorGMMStatistics`.
  * dropped exclusion of bipartite for sienaRI (only continuous excluded),
    but only if there are fewer second-mode nodes than actors.
    This required changes in
    `StatisticCalculator::calculateNetworkEvaluationStatistics`
    and in `siena07internals::getChangeContributionStatistics'.

# RSiena 1.3.7 (2022-02-18)

## Effects

  * New effect `avDeg`.

# RSiena 1.3.6 (2022-02-16)

## Effects

  * New effects  `simAllNear`,`simAllFar`, `absOutDiffIntn`, `avDegIntn`.
  * New effects `recipRateInv`, `recipRateLog` (Steffen Triebel).
  * Default internal effect parameter for `outOutActIntn`, `outOutAvIntn`,
    and `both` changed from 2 to 1.

## Functionality

   * Function `includeInteraction` now also can modify the `initialValue`
     of an effect; and the order of parameters for this function was changed,
     bringing it in line with `setEffect`.
   * Small clarifications of help pages for `includeInteraction` and
     `setEffect`.

# RSiena 1.3.5 and before

## Earlier changes

* See file `ONEWS_gf` in the source code at GitHub for changes
  in versions 1.2-26 to 1.3.5.

# RSiena 1.2-25 and before

## Earlier changes

* See file `ONEWS` in the source code at GitHub for changes
  in versions 17 to 1.2-25 (when the code was hosted at `R-forge`).
