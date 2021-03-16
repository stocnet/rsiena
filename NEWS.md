# RSiena 1.2.32

## 

2021-03-16

## Changes in RSiena:

### Effects:
   * New effects: homXTransRecTrip, toU.
   * This implied creation of a new effect class dyadANetNetObjective.
   * sqrt versions for parameter 2 for the effects to, toBack, toRecip,
     from, fromMutual.
   * Effects to, toU, toBack, toRecip, MixedInXW are dyadic.
   * Reinstated effect MixedInXW, also with sqrt version for parameter 2.
   * Dropped effect to.2 (identical to "to") 
     and MixedInWX (identical to "toBack").

### Improvements of functionality:
   * effectsDocumentation now also includes gmm effects (at the bottom).
   * Improved fromObjectToLaTeX in meta.table and siena.table.
   * Display of deviations from targets changed to after subtraction of targets.
   * Stop if no parameters are estimated and simOnly is FALSE (initializeFRAN).

### Reduction of functionality:
   * Vignette basicRSiena.Rmd dropped (available at website).

### Documentation:
   * Extended description of GMoM in the manual.
   * Description of toBack and toRecip in manual.
   * Changed keyword for some help pages.

### Corrections / safeguards
   * Correction in phase3.2 of a bug that sometimes led to an error message 
     if simOnly.
   * oneModeNet in effects.r: some further cases where the comparison of
     types with 'behavior' is replaced by 
     comparison with c('behavior', 'continuous').
   * Extra check in phase1.2.

### More neat code:
   * Dropped MixedOutStarFunction, MixedInStarFunction, MixedTwoPathFunction,
     (their functionality replaced by MixedTwoStepFunction).
   * Dropped MixedTwoStepFunction from effects 
     (its place is in effects\generic, and that's were it is).


# RSiena 1.2-31

## 

2021-02-27

## Changes in RSiena:
   * Generalized method of Moments implemented (Viviana Amati):
     see docs\manual\Changes_RSiena_GMoM.tex;
     new function includeGMoMStatistics, extended functionality of siena07.
   * Require R >= 3.5.0.
   * xtable added to "Imports" (used to be in "Suggests").
   * dyadicCov made to accept also changing dyadic covariates.
   * Used 'verbose' condition in sienaGOF also for last console output.
   * new arguments plotAboveThreshold and verbose for funnelPlot.

# RSiena 1.2.30

## 

2021-02-23

## Changes in RSiena:
   * Resolved issue with continuous dependent behavior variables 
     (Nynke Niezink).

# RSiena 1.2-29

## 

2020-12-10

## Changes in RSiena:
   * New effects (due to Christoph Stadtfeld):
     transtrip.FR, transtrip.FE, transtrip.EE, WWX.EE, WWX.FR, WXX.FE,
     WXX.ER, XWX.ER, XWX.FE, to.2, toBack, toRecip.
   * New effect transtripX.
   * New functions meta.table and funnelPlot.
   * For effect from.w.ind, option parameter=-1 added.
   * The to effect is an ego effect.
   * New parameter 'tested' in sienaGOF.
   * For siena.table, some of the effectNames changed to nice strings,
     so that LaTeX can run without errors if type='tex'.
   * The object produced by siena08 now has IWLS estimates more easily 
     accessible, as object$muhat and object$se.muhat.
   * Error message in sienaTimeTest for sienaFit objects produced with
     lessMem=TRUE.
   * More extensive error message for error in named vectors in algorithm object
     (checkNames in initializeFRAN).
   * For sienaDataCreate: more extensive error message, and class(...) replaced
     by class(...)[1]. 
   * multiplication factor added to print.sienaAlgorithm if maxlike.
   * In sienaAlgorithmCreate: requirements for mult corrected in help page.
   * In sienaAlgorithmCreate, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set. 


# RSiena 1.2-28

## 

2020-09-30

## Changes in RSiena:
   * Adapted filter "disjoint" so that it operates correctly
     also when the network is symmetric.
     Consequence: constraint that two networks are disjoint
     operates correctly also when one of the networks is symmetric
     and the other is not.
   * Adapted filter "higher" so that it operates correctly
     also when the other network is symmetric.
     Consequence: constraint that one network is at least as high
     as another network operates correctly also when 
     the higher networks is symmetric and the other is not.
   * In "CheckConstraints", used in "sienaDataCreate", the requirement 
     was dropped that the two networks have the same symmetry property;
     and for "higher" it is required that if the lower network
     is symmetric, the higher network is also symmetric.
   * In "sienaDataConstraint", if type is "disjoint" or "atLeastOne",
     the constraint is also implemented for the pair (net2, net1).
   * Vignette basicRSiena added (was earlier available as a script);
     thanks to James Hollway.

# RSiena 1.2-26

2020-09-17 

## Changes in RSiena:
   * Changed requirement for tcltk to "Suggests", 
     and further modified / cleaned up DESCRIPTION.
   * In siena07: if (!requireNamespace("tcltk")) set batch to TRUE.
   * In NAMESPACE drop tcltk
   * In sienaAlgorithmCreate, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set.


