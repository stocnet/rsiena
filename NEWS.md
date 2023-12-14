# RSiena 1.4.2

2023-12-14

## Changes in RSiena:

### Effects:
  * New effect `outXMore`, `outMore3`. 
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
   
# RSiena 1.3.5

##

2021-12-15

## Changes in RSiena:  

### Bug corrections  
   * Corrected the check for effects in `initializeFRAN.r` which led to errors
     if interaction effects are included, because of the changes to
     `includeInteraction` in version 1.3.4.

# RSiena 1.3.4

##

2021-12-08

## Changes in RSiena:  

### Effects:
  * New effects `inRateInv`, `inRateLog` (Steffen Triebel).

### Improvements of functionality:  
   * When an effects object with interaction effects is printed,
     the names of the interacting effects are mentioned,
     and prefixes "int." and "i3." were dropped.
   * The check of whether an interaction effect is allowed now is done
     immediately when creating the interaction effect instead of waiting
     for its use in `siena07`. 
   * For function `sienaGroupCreate` some changes were made:   
     if it is applied to a list of length 1, attributes of the single group
     are not recomputed;   
     if it is applied to a list of length larger than 1, the attributes
     "range" and "range2" of behavioral variables of individual groups are
     computed as the range of the unions of the ranges of all the groups.
     The same is done for covariates.
   * `print.sienaGroup` slightly extended.
   * Creation of covariates gives a warning (optional) if all values 
     are missing, and also if all non-missing values are the same.

### Bug corrections  
   * if a `sienaGroup` object is given to `sienaBayes` and some of the 
     covariates are constant in one or more of the groups, the ``simX` effect 
     will not run into an error any more; this is achieved by the first change
     mentioned above for `sienaGroupCreate`.
   * `sienaFitThetaTable` in `sienaprint.r` was corrected for `from sienaBayes`.
   * `siena.table` was corrected for `sienaBayes` objects, and this possibility
     was mentioned in the help file.

### Coding changes:
   * In `sienaprint.r`, methods and functions relating to `sienaBayes` omitted.

### Tests:   
   * Test 18, a test for ``sienaGroupCreate`, added to `parallel.R`.
     

# RSiena 1.3.3
   
## 
   
2021-10-09

## Changes in RSiena:  

### Effects:
   * internal effect parameter for diffusion rate effects ('at least p').
   * New effect `outOutAvIntn`.
   * `outOutActIntn` also made available for non-directed explanatory networks 
     and for two-mode dependent networks.

### Improvements of functionality: 
   * `toggleProbabilities` added to output of `sienaRI`.
   * Trial values of `theta` used during Phase 2 of `siena07` added 
     to the `sienaFit` object `ans` as `ans$thetas`.
   * warnings if a data object contains only missings, or only the same value.

### Bug corrections  
   * if a data set contains a constant covariate, the simX effect will
     not run into an error any more; this is obtained by defining the `range` 
     attribute of the covariate as 1 to prevent division by zero,
     and the `simMean` attribute as 0 instead of `NaN` (`sienaDataCreate`). 
     This is relevant especially for `sienaGroup` data sets, where covariates 
     might be constant for some of the groups.`

### Corrections of documentation: 
   * In `sienaAlgorithmCreate`, the default values of `diagonalize` for MoM
     is 0.2; this was corrected in the help file.
     

# RSiena 1.3.2
   
## 
   
2021-07-29

## Changes in RSiena:  

### Improvements of functionality: 
   * Effects of type `creation` or `endow` represented in siena.table
     by `creation` and `maintenance`, respectively.    
     This was erroneously omitted in 1.3.1.

# RSiena 1.3.1
   
## 
   
2021-07-27

### Effects:
   * New effects: `crprodInActIntn` (Nynke Niezink), `XXW`.

### Improvements of functionality: 
   * Effects of type `creation` or `endow` represented in siena.table
     by `creation` and `maintenance`, respectively.
   * `updateTheta` also accepts `sienaBayesFit` objects as `prevAns`.

## Small corrections:    
   * If `upOnly` or `downOnly`, the (out)degree (density) effect is also 
     excluded for symmetric networks 
     (this was reported in `print01Report`, but not carried out).
     This happens in `effects.r`.
   * Message corrected in `sienaDataCreate` if there is an attribute `higher`.
   
   
# RSiena 1.3.0
      
## 
      
2021-05-02  

   * Drop `testthat` in tests.   

# RSiena 1.2.34
   
## 
   
2021-04-30  
   
### New functions:
   * `testSame.RSiena`.
   
### Effects:  
   * New effects: `avInSim` (thanks to Steffen Triebel), `totInSim`, 
     `avInSimPopAlt`, `totInSimPopAlt`, `constant`,
     `avAttHigher`, `avAttLower`, `totAttHigher`, `totAttLower`.
   * Changed effects: endowment and creation types for `avInSim`
     (brought in line with these types for `avSim`).
    
### Improvements of functionality:  
   * `funnelPlot` adapted to lists of `sienaFit` objects
     containing missing estimates or standard errors.
   * `plot.sienaGOF`: new parameter `position`.
   * Small improvements (length of effect names) in `meta.table` and 
    `siena.table`.  
   
### Bug corrections  
   * Restore backward compatibility with respect to checks of `x$gmm`.  
   * In test functions: correct names reported for tested effects by 
     using `ans$requestedEffects` instead of `ans$effects`.
   
### Code improvements   
   * Improved coding of `SimilarityEffect`, using new parts
     of `NetworkDependentBehaviorEffect`.   
   * Changed unsigned actors to int in `Continuousvariable` and 
     `EpochSimulation`;   
     int ...EffectCounts to unsigned in `BehaviorVariable`, 
     to avoid warnings in C++ compilation.
   * Changed name of `similarity(int i, int j)` to `actor_similarity`
     in order to avoid confusion with `similarity(double v1, double v2)`.
   
### Corrections  
   * Took out of `NAMESPACE` a few imported functions from `graphics`, 
     `stats`, `utils` that were not used.  
   * Correction of footer of `CovariateDistance2EgoAltSimNetworkFunction.h`.
  
# RSiena 1.2.33

## 

2021-03-19  


   * Adjusted `configure`, `cleanup` and `Makevars` files for just C++ checks.
   * Pandoc dropped as a system requirement.
 
# RSiena 1.2.32

## 

2021-03-16

### Effects:
   * New effects: `homXTransRecTrip`, `toU`.
   * This implied creation of a new effect class `dyadANetNetObjective`.
   * sqrt versions for parameter 2 for the effects `to`, `toBack`, `toRecip`,
     `from`, `fromMutual`.
   * Effects `to`, `toU`, `toBack`, `toRecip`, `MixedInXW` are dyadic.
   * Reinstated effect `MixedInXW`, also with sqrt version for parameter 2.
   * Dropped effect `to.2` (identical to `to`) 
     and `MixedInWX` (identical to `toBack`).

### Improvements of functionality:
   * `effectsDocumentation` now also includes `gmm` effects (at the bottom).
   * Improved `fromObjectToLaTeX` in `meta.table` and `siena.table`.
   * Display of deviations from targets changed to after subtraction of targets.
   * Stop if no parameters are estimated and `simOnly` is FALSE (initializeFRAN).

### Reduction of functionality:
   * Vignette `basicRSiena.Rmd` dropped (available at website).

### Documentation:
   * Extended description of GMoM in the manual.
   * Description of `toBack` and `toRecip` in manual.
   * Changed keyword for some help pages.

### Corrections / safeguards
   * Correction in phase3.2 of a bug that sometimes led to an error message 
     if `simOnly`.
   * `oneModeNet` in `effects.r`: some further cases where the comparison of
     types with `behavior` is replaced by 
     comparison with c(`behavior`, `continuous`).
   * Extra check in `phase1.2`.
   * Temporarily drop the final part of `test16`,
     in view of an irreproducible error.

### More neat code:
   * Dropped `MixedOutStarFunction`, `MixedInStarFunction`, 
     `MixedTwoPathFunction`, 
     (their functionality replaced by `MixedTwoStepFunction`).
   * Dropped `MixedTwoStepFunction` from effects 
     (its place is in `effects\generic`, and that's were it is).


# RSiena 1.2-31

## 

2021-02-27

   * Generalized method of Moments implemented (Viviana Amati):
     see docs\manual\Changes\_RSiena\_GMoM.tex;
     new function `includeGMoMStatistics`, extended functionality of `siena07`.
   * Require R >= 3.5.0.
   * `xtable` added to `Imports` (used to be in `Suggests`).
   * `dyadicCov` made to accept also changing dyadic covariates.
   * Used `verbose` condition in `sienaGOF` also for last console output.
   * new arguments `plotAboveThreshold` and `verbose` for `funnelPlot`.

# RSiena 1.2.30

## 

2021-02-23

   * Resolved issue with continuous dependent behavior variables 
     (Nynke Niezink).

# RSiena 1.2-29

## 

2020-12-10

   * New effects (due to Christoph Stadtfeld):
     `transtrip.FR`, `transtrip.FE`, `transtrip.EE`, `WWX.EE`, `WWX.FR`, 
     `WXX.FE`, `WXX.ER`, `XWX.ER`, `XWX.FE`, `to.2`, `toBack`, `toRecip`.
   * New effect `transtripX`.
   * New functions `meta.table` and `funnelPlot`.
   * For effect `from.w.ind`, option parameter=-1 added.
   * The `to` effect is an ego effect.
   * New parameter `tested` in sienaGOF.
   * For `siena.table`, some of the effectNames changed to nice strings,
     so that `LaTeX` can run without errors if `type=tex`.
   * The object produced by `siena08` now has IWLS estimates more easily 
     accessible, as `object$muhat` and `object$se.muhat`.
   * Error message in `sienaTimeTest` for `sienaFit` objects produced with
     `lessMem=TRUE`.
   * More extensive error message for error in named vectors in algorithm object
     (`checkNames` in `initializeFRAN`).
   * For `sienaDataCreate`: more extensive error message, and `class(...)` 
     replaced by `class(...)[1]`. 
   * multiplication factor added to `print.sienaAlgorithm` if `maxlike`.
   * In `sienaAlgorithmCreate`: requirements for `mult` corrected in help page.
   * In `sienaAlgorithmCreate`, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set. 


# RSiena 1.2-28

## 

2020-09-30

   * Adapted filter `disjoint` so that it operates correctly
     also when the network is symmetric.
     Consequence: constraint that two networks are disjoint
     operates correctly also when one of the networks is symmetric
     and the other is not.
   * Adapted filter `higher` so that it operates correctly
     also when the other network is symmetric.
     Consequence: constraint that one network is at least as high
     as another network operates correctly also when 
     the higher networks is symmetric and the other is not.
   * In `CheckConstraints`, used in `sienaDataCreate`, the requirement 
     was dropped that the two networks have the same symmetry property;
     and for `higher` it is required that if the lower network
     is symmetric, the higher network is also symmetric.
   * In `sienaDataConstraint`, if type is `disjoint` or `atLeastOne`,
     the constraint is also implemented for the pair (net2, net1).
   * Vignette `basicRSiena` added (was earlier available as a script);
     thanks to James Hollway.

# RSiena 1.2-26

2020-09-17 

   * Changed requirement for `tcltk` to `Suggests`, 
     and further modified / cleaned up DESCRIPTION.
   * In siena07: if `(!requireNamespace(tcltk))` set batch to TRUE.
   * In NAMESPACE drop `tcltk`.
   * In `sienaAlgorithmCreate`, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set.


