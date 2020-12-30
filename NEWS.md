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
   * SienaAlgorithmCreate: requirements for mult corrected in help page.
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


