2025-10-28

# RSiena 1.5.6

## Changes in RSiena:
### Coding
  * The virtual function `preprocessEgo` now is defined properly for 
    the function classes `NetworkAlterFunction`, 
   `CovariateNetworkAlterFunction`, `DoubleCovariateFunction`, 
   `SameCovariateInTiesFunction`. 
   This was done to correct a bug that appeared for effect `sameXInPop`
   for a two-mode network; perhaps there was a similar bug for `diffXInPop`,
   `sameXOutAct`, `diffXOutAct`, and `sameXVInPop`. 
### Functionality:
  * Auxiliary function `egoAlterCombi` for GoF now omits values 
    where the behavioral variable is `NA`.
  * Various small additions to output of `siena07` for `verbose=TRUE` 
    (`phase1.r` and `phase2.r`). 

2025-10-17

# RSiena 1.5.3

### Effects
  * New effects `quad_cc`, `avAlt_cc`, `totAlt_cc`.

2025-09-06

# RSiena 1.5.2

## Changes in RSiena:
### Functionality:
  * `selectionTable` got an attribute `quad` indicating whether the plot
    is a quadratic function.
  * Users can now extract now extract the changeContributions when running
    `siena07` by setting an argument `returnChangeContributions=TRUE`. 
    If used together with `nsub=0` and `prevAns` or modified initial 
    values in the effects object, especially useful for post-estimation,
    e.g. in `sienaRIDynamics`.
  * `sienaRIDynamics` uses siena07 directly now and is reinstated.
### Effects
  * New effects `outActMore_ego`, `outActSqrtMore_ego`, `outMore_ego`,
    `outPopMore`, `outPopSqrtMore`, `outPopThreshold`.
### Coding
  * Use `(any(!gmm))` in  `initializeFRAN` to allow the use of a `prevAns` 
    object with a different estimation method for a multigroup estimation.

2025-07-12

# RSiena 1.5.1

## Changes in RSiena:
### Functionality:
  * `sienaRI` reinstated.
  * Option `prML=2` for maximum likelihood estimation 
    using the `move` proposal step reinstated (`sienaAlgorithmCreate`).
### Effects
  * New effects `divOut_ego`, `divIn_ego`. 
### Coding
  * Imported functions from packages `Matrix`, `lattice`, `parallel`, `MASS`, 
    and `xtable` mentioned specifically in the `Namespace`
    instead of importing these entire packages.


# RSiena 1.5.0

2025-07-05

## New CRAN version

## Changes in RSiena:
### Maintainer
  * Christian Steglich now is maintainer.

# RSiena 1.4.25

2025-07-05

## Changes in RSiena:
### Coding
  * Undid other changes in version 1.4.23 to
    `getTargetsChangeContributions` in `siena07setup.cpp`.

# RSiena 1.4.24

2025-07-04

## Changes in RSiena:
### Effects
  * Error for creating effects in effects group `dyadSecondBipartiteObjective` 
    was corrected in file `effects.r`, 
    so they now are included also for changing dyadic covariates.
    This affected effects `XWX`, `XWX1`, and `XWX2`.
### Coding
  * Function `Chains:printConsecutiveCancelingPairs` deleted from `model\ml`,
    because it was not used and led to a protection error.
  * Undid changes in version 1.4.23 to `getTargetActorStatistics` and 
    `getTargetsChangeContributions` in `siena07setup.cpp`.

# RSiena 1.4.23

2025-05-03

## Changes in RSiena:
### Bug corrections
  * Repaired (hopefully) memory leak in `getTargetActorStatistics` and 
    `getTargetsChangeContributions` in `siena07setup.cpp`. 
  * Repaired (hopefully) memory leak in `move` in `MLSimulation.cpp`. 
  * Repaired memory leak  in `DoubleCovariateCatFunction`. 
### Effects
  * New effects `fromAny`, `sameInXCycle4`. 
  * Effect shortName `cycle4ND` replaced by `cycle4`.
  * Better treatment of missing covariate values in effect `sameXCycle4`. 
  * Internal effect parameter values 3 and 4 for `sameXInPop`,`diffXInPop`,
    `sameXInPopIntn`,  `sameXInActIntn`, `homXOutAct2`. 
  * New C++ class `CatCovariateDependentNetworkEffect` (for `homXOutAct2`). 
  * Added `(#)` to the `effectName` of `from`, `sameXInPop`, `diffXInPop`, 
	`sameXVInPop`, `sameXVInPop2`.
### New functionality
  * In `print.sienaEffects(..., includeShortNames=TRUE)`, the `effectNumber`
    is also printed.
### Coding
  * Class `Covariate` has new variable `covariateN`, which then is transferred
    as `covarN` to `CovariateDependentNetworkEffect` and 
    `DoubleCovariateFunction`. 
### Changes in documentation:
  * Definition of effects `outRateLog`, `inRateLog`, and `recipRateLog`
    corrected in the manual.

# RSiena 1.4.22

2025-02-08

## Changes in RSiena:
### New functionality
  * New auxiliary GOF function `egoAlterCombi`. 
  * Parameter `showAll` added to `plot.sienaGOF`.
### New src functionality
  * New table `IntLogTable` and new generic function`IntLogFunction`. 
### Effects
  * Internal parameter 0 (for log(x)) added for `outActIntn`.
    (It would be trivial to implement this also for the other
    mixed degree effects, but currently there seems no need.)
### Bug correction.
  * In `sienaGOF`,  if the auxiliaryFunction does not always 
    give vectors of the same length, the error message gives properly
    the name of the auxiliaryFunction.

# RSiena 1.4.21

2024-12-18

## Changes in RSiena:
### New functionality
  * New functions `selectionTable` and `influenceTable`, 
    with `print` methods. 
### Effects
  * New effect `altHigherEgoX`.
  * Effect `higher` also implemented for symmetric networks.
  * Used # in name of threshold effects.
### Bug corrections
  * In `siena.table`,  fixed parameter values are not reported as NA,
    but as their fixed values.
### Messages
  * If `includeInteraction` is called with argument `parameter`, 
    an error message appears that this keyword should not be given. 
  * Extension of message in `sienaTimeTest` given in the case of
    collinearities. 
### Help pages
  * Improved explanation of the use of `interaction1` and `interaction2`
    in the help page for `includeInteraction`. 

# RSiena 1.4.20

2024-11-10

## Changes in RSiena:
### Effects
  * New effects `varAlt` and `avSimVarAlt`.
### Bug corrections
  * Interaction effects for continuous behavior now turned off. The R side of
    `includeInteraction` for continuous behavior had been implemented but 
    not the C++ side, so now an error message is given.
### Error messages
  * If `includeInteraction' is called for a continuous behavior effect, an
    error message is given that interaction effects are not yet implemented
    for continuous behavior.
### Coding
  * For the new effects, `lvariance`, the behavior variance over all but the 
    last wave, was included as member data of `BehaviorLongitudinalData` 
    together with the corresponding member functions. 

# RSiena 1.4.19

2024-09-03

## Changes in RSiena:

### Effects
  * New effects `outThreshold` and `outThreshold2`.
  * `outActIntn` is an ego effect (`allEffects.csv`).
### New functionality:
  * Network option 11 = `NETCONTEMP` for use of contemporaneous statistics
    for estimating all evaluation effects of the network variable.
### Bug corrections
  * Error corrected that occurred for estimation by the GMoM if
    some effects are fixed.
  * List of effect names when `thetaBound` is exceeded corrected (`phase2.r`).
### Error messages
  * If `setEffect` is called for an effect with `type=gmm`, an error message
    is given that the function called should be `includeGMoMStatistics`. 

# RSiena 1.4.18

2024-07-31

## Changes in RSiena:

### Effects
  * New gmm effects for co-evolution of multiple networks: 
    `crprod_gmm`, `to_gmm`, `from_gmm`.
### Bug corrections
  * Bug corrected in `sienaGOF` which occurred for models with tested effects
    if `iterations` is less than `sienaFitObject$n3`. 
### Coding
  * The new gmm effects required the extension of the function classes
    `AlterFunction`, `NetworkAlterFunction`, and `MixedNetworkAlterFunction`
    to allow making estimation statistics depend totally on the simulated state.
  * Function `fixUpEffectNames(effects, defaultEffects)`
    moved from `initializeFRAN.r` to `sienaEffects.r`. 

# RSiena 1.4.13

2024-06-06

## Changes in RSiena:

### Effects:
  * New effect `homXOutAct2`.
### Bug corrections: 
  * Function `setEffect` corrected (it did not give the proper internal
    parameters in the `effectName`). 
  * Function `updateSpecification` corrected (it did not work properly
    for including interactions).
### Improved functionality:
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


# RSiena 1.4.12

2024-04-27


## Changes in RSiena:

### Improved coding:
  * All objects created by functions, if not print or summary,
   now have an attribute "version", which is the package version.
### Miscellanea:
  * New RSiena logo inst/rsienalogo-2.png


# RSiena 1.4.11

2024-04-23


## Changes in RSiena:

### Effects:
  * Effect name `outOutDist2AvIntn` changed to `avAlt.2M.tot`.
  * New effects `avAlt.2M.tie`, `avAlt.2M.tot`, `avAltU.2M.tie`, 
    `dist2OutInActIntn`, `nDist2ActIntn`, `sharedToU`.
  * Changed sqrt treatment in `outOutDist2ActIntn` and
    `outOutDist2AvIntn` / `avAlt.2M.tot`.
### New functionality:
  * Method `print.sienaEffects` has an extra parameter `includeShortNames`
    to do what the name of this parameter suggests.
  * Improved error message for function `updateSpecification`.

# RSiena 1.4.10

2024-03-25


## Changes in RSiena:

### Effects:
  * Effects `sameXV` and `sameXVInPop` added for symmetric networks, 
    and restricted to integer-valued variables in the range from 0 to 20.
### New functionality:
  * Parameter `silent` (new in version 1.4.8) 
    in `sienaAlgorithmCreate` activated. 
### Improved coding:
  * Actor covariates in `sienaData` have a new attribute `lowIntegers` 
    used for in/excluding effects `sameXVInPop` and `sameXV` in
    `getEffects`. 

# RSiena 1.4.9

2024-03-21


## Changes in RSiena:

### Effects:
  * New effects `crossXOutAct`, `outOutDist2ActIntn`, 
    `outOutDist2AvIntn`, `inPopOutW`.
  * New effect group `doubleCovarNetObjective`.
  * New effects `sameXV` and `sameXVInPop` for bipartite networks.
  * `sameXCycle4` added for one-mode and symmetric networks.
  * `sharedTo` gets default internal effect parameter `p=3`. 
### Improved functionality:
  * Function `updateSpecification` now also updates interaction effects
    and `initialValues`.
### Improved coding:
  * Internal functions `numberIntn`, `numberBehIntn`, `checkVersion` 
    defined in file `initializeFRAN.r`.
  * The `Covariate` class and its descendants (all actor covariates)
    now have functions `min` and `max`. 

# RSiena 1.4.8

2024-02-29


## Changes in RSiena:

### Bug corrections:
  * Correction of memory leak in `siena07setup.ccp` for ML estimation.
### New functionality:
  * New parameter `silent` in `sienaAlgorithmCreate`. 

# RSiena 1.4.7

2024-02-21

## New CRAN version

## Changes in RSiena:

### Bug corrections:
  * Correction (again) of error message in `siena07utilities::Rterminate`.

# RSiena 1.4.6

2024-02-19


## Changes in RSiena:

### New functionality:
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
### Bug corrections:
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
### Dropped functionality:
  * `sienaRI` temporarily disabled because of a memory leak.
  * The option `prML=2` temporarily disabled because of a memory leak in the
    `move` proposal distribution (`sienaAlgorithmCreate`).
  * `doMoreUpdates` moved from `maxlike.r` to `maxlikecalc.r`.
    The rest of `maxlike.r` as wel as `maxlikefn.Rd` deleted.
    These were not used anywhere.
### Other changes:
  * List of changes in the code before 2022 moved from `NEWS.md` to `ONEWS_gh`.

# RSiena 1.4.5

2024-02-14

## Changes in RSiena:

### Package

  * Migrated package repository to "stocnet" organisation.

### Bug corrections:
  * Put `#include <Rinternals.h>` as the last of the include commands
    in various `.cpp` files (Tomas Kalibera). 

# RSiena 1.4.2

2023-12-14

## Changes in RSiena:

### Effects:
  * New effects `outXMore`, `outMore3`. 
  * `Interactiontype` of `altLThresholdX` and `altRThresholdX` is dyadic.
  * `Interactiontype` of `degAbsDiffX`, `degPosDiffX`, and `degNegDiffX` 
    is "" (blank). 
  * Corrected effect `outMore`. 
### Improved coding:
  * `bxeffects` initialized to 0 in `ContinuousVariable::accumulateScores`.
  * All private variable declarations in the C++ `*.h` files
    were initialized using braces. 
  * In `mixedTriadCensus`, the check of the condition about the matrix 
    dimensions was split into its two parts.
  * `TruncatedOutdegreeEffect2` dropped from `src/model/effects` 
    (was superfluous).
### New functionality:
  * New parameter `iterations` in `sienaGOF` to allow shorter computations.

# RSiena 1.4.1

2023-11-01

## New CRAN version.

### Changes in help files:
  * Changes in accordance with "Guidelines for Rd files". 

# RSiena 1.4.0

2023-10-31

## Changes in RSiena:

### Changes in meta-data:
  * `Author` field omitted in `DESCRIPTION`, because `Author@R` is sufficient.
### Changes in `inst` directory:
  * Superfluous files in `inst` directory omitted.
  * New `CITATION`.
### Small changes in coding:
  * Superfluous "lsimulatedDistance" in `SdeSimulation.h` deleted.
  * In `PrimarySettingEffect.cpp`, used `to_string` for converting a number
    to string in an error message.

# RSiena 1.3.28

## 
   
2023-10-11

## Changes in RSiena:

### Changed effect:
  * `threshold`, `threshold2`, `threshold3`, `threshold4` changed to 
    work with non-centered parameters (not backward-compatible).
### Improved documentation:
  * Descriptions of effects `altInDist2W` and `totInDist2W` added
    to the manual (the effects had been there since a long time, 
    but not documented). 

# RSiena 1.3.27

## 
   
2023-09-29

## Changes in RSiena:

### Coding:
  * Corrected one line in `siena07models.cpp`, which led to slowness of 
    `siena07` since version 1.3.18.
### New effects:
  * `avInSimDist2`, `totInSimDist2`, `sameEgoDist2`,`sameEgoInDist2`,
    `outMore2`, `divOutEgoIntn`, `divInEgoIntn`, `divOutAltIntn`, 
    `divInAltIntn`.
  * `avTAltDist2` and `totTAltDist2` also implemented for behavior co-evolving
    with symmetric networks.
### Documentation:
  * Some explanation is given in the manual about internal effect parameters
    for interactions created by `includeInteraction`.

# RSiena 1.3.26

## 
   
2023-08-15

## Changes in RSiena:

### Coding:
  * Improved Phase 1 derivative matrix computation for basic SDE parameters.
  * Added continuous behavior to returned simulated data.
### Corrections:
  * Period/groupwise tests in `sienaTimeTest` corrected for the case of 
    non-saturated sets of dummy variables.
  * `plot.sienaTimeTest` for "pairwise=TRUE" changed so that the warning
    is avoided. 
  * `sienaGOF` corrected so that again it can handle auxiliary functions
    referring to more than one `varName` (such as in `mixedTriadCensus`).


# RSiena 1.3.24

## 
   
2023-08-01

## Changes in RSiena:

### Corrections:
  * In `getEffects`, the effects object was given an attribute `version`,
    which was not done correctly in version 1.3.23.
    (This led to always giving a warning if any interaction effects were 
    specified.)
  * Corrections of implementation of acceptance by `sienaGOF` of a list of 
    `sienaFit` objects (was not correct in version 1.3.23). 
### Additional testing:
  * function `includeInteraction` used in "parallel.R" (for testing). 


# RSiena 1.3.23
   
2023-06-29

## Changes in RSiena:

### New effects: 
  * New effects `diffWXClosure`, `sameWWClosure`,  `diffWWClosure`, 
    `diffXWClosure`, `sameXWClosure`, `unequalX`.
  * `JoutMix` made available for bipartite dependent networks.
  * For continuous behavior variables depending on a bipartite dependent 
    network, the effect group `continuousBipartiteObjective` was created,
    with effects `outdeg`, `outdegSqrt`, and `isolateOut`. 
  * `sameXOutAct` and `diffXOutAct` now have a parameter 2 for `sqrt`.
### Corrections:
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
### New functionality:
  * The model for continuous behavior variables seems to work now,
    because of the first correction mentioned above.
  * `sienaGOF` now also accepts a list of `sienaFit` objects.
### Improved coding:
  * Better text for stop in `initializeFRAN` when there is a mismatch
    between effects objects disabling the creation of interaction effects.
  * Warning in `initializeFRAN` if the version of the effects object 
    is not current and the effects object contains interaction effects
    (then it is possible that the interacting effects are chosen incorrectly,
    even though the `effectName` of the interaction seems OK).
  * Better error message in `sienaGOF` if `groupName` or `varName` is incorrect.
  * Use default bandwidth selection in violin plot for `sienaGOF`
    (the use of "nrd" sometimes led to absent plots because of negative bw). 

# RSiena 1.3.22

## 
   
2023-05-11

## Changes in RSiena: 
### Coding:
  * Corrected and cleaned up virtual definitions in `AlterFunction` and its
    descendants, in particular `CovariateNetworkAlterFunction`.
  * Added `const` to virtual specification of `value` in `AlterFunction.h`
    and all of its descendants.
  * Replaced ambiguous call to `std::abs` in `AbsDiffFunction.cpp`. 
### New functionality:
  * For one-mode networks, new model options `DOUBLESTEP25`, `DOUBLESTEP50`,
    `DOUBLESTEP75`, `DOUBLESTEP100`. 
### Corrections:
  * The first item in "Coding" implies correction of several distance-2 network 
	 effects such as `altDist2`, `totDist2` and `altInDist2`.
  * In `sienaAlgorithmCreate`, changed default `prML=2` back to `prML=1`; 
    stop if Maximum Likelihood estimation is attempted for a data set
    containing more than one dependent variable 
    with `prML=2` (implemented in `initializeFRAN.r`). 

# RSiena 1.3.20

## 
   
2023-04-22

## Changes in RSiena:  
### Corrections:
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
### Improvements of functionality:
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
### Documentation:
  * Reference about score-type test added to `Wald.Rd`.
  * In the help page for `sienaDependent`, it is mentioned that if there are 
    one-mode as well as two-mode dependent networks,
    the one-mode networks should come first.

# RSiena 1.3.19

## 
   
2023-02-07

## Changes in RSiena:  
### Coding:
  * `siena07internals.cpp` adapted to be compatible with new clang 16 C++ 
    compiler (thanks to Brian Ripley).
### New effects:
  * New effect `inPop_dya`.
  * Parameter 2 for `sameXInPop` and `diffXInPop`.
### Corrections:
  * Help page for `siena07` corrected with respect to `x$lessMem`.
### Improvements of functionality:
  * `coCovar` and `varCovar` now can handle variables with only one 
    non-missing value, but will stop with an error message 
    if all values are missing.

# RSiena 1.3.18

## 
   
2023-01-29

## Changes in RSiena:  
### Improvements of functionality:
  * Additional step type `move` for MH proposal distribution
    for likelihood estimation (thanks to Charlotte Greenan).
  * Accordingly, parameters changed that are used in `sienaAlgorithmCreate`
   for probabilities of MH steps, now summarized in `prML`; with a new default.
  * List elements `accepts`, `rejects`, `aborts` for `sienaFit` objects
    produced by ML estimation improved/corrected by reorganizing them in C++.
  * List element `ac3` added to `sienaFit` object if `maxlike`.

# RSiena 1.3.17

## 
   
2023-01-06

## Changes in RSiena:  

### Improvements of functionality:
  * `sienaGOF` now accepts simulated auxiliary statistics containing missing
    values. If there are any, this will be reported with a warning
    if `giveNAWarning` is `TRUE`.
  * `sienaDataCreate` now also accepts, as "...", a list of such objects.

# RSiena 1.3.16

## 
   
2023-01-02

## Changes in RSiena:  

### Corrections:

### Effects:
    `inPopIntnX`, `inActIntnX`, `outPopIntnX`, `outActIntnX`, `sameXInPopIntn`, 
    `sameXOutPopIntn`, `sameXInActIntn`, `sameXOutActIntn` restored
    (these had got lost in some way...).
### Updates:
  * All occurrences of `http` in `R` and `Rd` files changed to `https`.
  * `seq_len` used and superfluous `c()` omitted in various R files.

# RSiena 1.3.15

## 
   
2022-11-27

## Changes in RSiena:  

### Corrections:
  * `siena08`: correct p-value `pTsq` for overall test statistic `Tsq`
  * `print.summary.sienaMeta`, `siena07`, `print01Report`: drop RForge revision.
  * Correct "objname" to "projname" in `meta.table` (`siena08.r`).
  * Simplify `LaTeX` output of `meta.table`. 
  * `seq_along` and `seq_len` used in `print01Report`. 

# RSiena 1.3.14

## 
   
2022-11-04

## Changes in RSiena:  

### Note:
  * CRAN version.

### Corrections:
  * Update `configure` and  `configure.ac` (with help from Brian Ripley). 

# RSiena 1.3.13

## 
   
2022-10-07

## Changes in RSiena:  

### Updates:
  * Replacements in EffectFactory.cpp of single | operator by ||.

# RSiena 1.3.12

## 
   
2022-10-06

## Changes in RSiena:  

### Updates:
  * Changes to comply with new version of `Matrix` package.
  * Replacements in some C++ functions of single & and | operators by && and ||.
### Corrections:  
  * `universalOffset` initialized as 0; it was earlier initialized as
    the maximum real number (`NetworkLongitudinalData.cpp`). 
  * `thetaStore` deleted (was trash in `phase2.r`).
  * Various comparisons for vectors with 0 changed to using `all`
    to avoid warnings (`initializeFRAN.r`).
### Code modifications:
  * `sigmas` and `meansigmas` added to `sienaRI` object.
  * Print of standard deviations in the `sienaRI` object for `printSigma=TRUE` 
    changed to using averages at the variance level.
  * If `returnThetas` in the call of `siena07`, also simulated estimation 
    statistics during Phase 2 (deviations from targets) are returned.
### Effects:
  * Several new effects related to primary setting:
    `nonPCompress`, `primCompress`, `primary`, `primDegAct`,
    `primDegActDiff`, `primDegActDiffSqrt`, `primDegActSqrt`,
    `primDegActLog`, `primDegActInv`.
  * `gwdspFB` effect added for two-mode networks.
  * New effects `outAct_ego`, `inAct_ego`,`reciAct_ego`, `toAny`.
  * For effects `to`, `toBack`, `toRecip`, `mixedInXW`, 
    internal effect parameter 3 now specifies truncation of the number of 
    twosteps (change to `MixedTwoStepFunction`). 
### Improvements of documentation:
  * Modified help page for `sienaRI`.
  * Small modifications of help page for `sienaGOF`.

# RSiena 1.3.11

## 
   
2022-05-30

## Changes in RSiena:  

### Corrections:   
  * Correction in `effects.r` of error that led to warning
    for multivariate networks.
  * Correction of help page for `sienaGOF` (`groupName`).
  * Correction of `igraphNetworkExtraction` in the help page for
    `sienaGOF-auxiliary`.

### Improvements of functionality: 
  * Further explanation of `mixedTriadCensus` in the help page for
    `sienaGOF-auxiliary`.


# RSiena 1.3.10

## 
   
2022-04-28

## Changes in RSiena:  

### Corrections:
  * Bug corrected that occurred when several two-mode networks were included
    in the dependent variables, with an order restriction between them.
    (Correction of `HigherFilter` and `DisjointFilter`).

### Effects:
  * New effects `avInAltW`, `avWInAlt`, `totInAltW`, `totWInAlt` 
    (with help from Robert Krause).
  * Corrected implementation of `sharedTo`.
  
### Code modifications:
  * Several modifications to enable traceback of errors occurring
    in checkSenderRange called in inTies.

# RSiena 1.3.9

## 
   
2022-03-18

## Changes in RSiena:  

### Effects:
  * Corrected implementation of `simAllNear` and `simAllFar`.

### Corrections:
  * small correction of `summary.sienaGOF`.
  * small correction of `sienaTimeTest`.

# RSiena 1.3.8

## 
   
2022-03-07

## Changes in RSiena:  

### Effects:
  * Changed default internal effect parameter for `simAllNear` to 2 and for
    `simAllFar` to 4.

### Improvements of functionality: 
  * In `sienaTimeTest`, added `warn=FALSE` to `varCovar()` to avoid warnings.
  * Small improvements to help pages for `sienaGroupCreate` and `sienaGOF`.

### Corrections:
  * corrected sienaRI for behavioral variables.
    This required changes in 
    `StatisticCalculator::calculateBehaviorStatistics` and
    `StatisticCalculator::calculateBehaviorGMMStatistics`.
  * dropped exclusion of bipartite for sienaRI (only continuous excluded),
    but only if there are fewer second-mode nodes than actors.
    This required changes in 
    `StatisticCalculator::calculateNetworkEvaluationStatistics`
    and in `siena07internals::getChangeContributionStatistics'.


# RSiena 1.3.7

## 
   
2022-02-18

## Changes in RSiena:  

### Effects:
  * New effect `avDeg`.

# RSiena 1.3.6

## 
   
2022-02-16

## Changes in RSiena:  

### Effects:
  * New effects  `simAllNear`,`simAllFar`, `absOutDiffIntn`, `avDegIntn`.
  * New effects `recipRateInv`, `recipRateLog` (Steffen Triebel).
  * Default internal effect parameter for `outOutActIntn`, `outOutAvIntn`, 
    and `both` changed from 2 to 1.

### Improvements of functionality:  
   * Function `includeInteraction` now also can modify the `initialValue`
     of an effect; and the order of parameters for this function was changed,
     bringing it in line with `setEffect`.
   * Small clarifications of help pages for `includeInteraction` and
     `setEffect`.

# RSiena 1.3.5 and before

## Earlier changes:

###
* See file `ONEWS_gf` in the source code at GitHub for changes 
  in versions 1.2-26 to 1.3.5.


# RSiena 1.2-25 and before

## Earlier changes:

###
* See file `ONEWS` in the source code at GitHub for changes
  in versions 17 to 1.2-25 (when the code was hosted at `R-forge`). 
   
