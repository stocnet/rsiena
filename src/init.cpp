
#include <stdlib.h> // for NULL
#include "siena07setup.h"
#include "siena07models.h"

// R stuff last redefine length macro
#include <R.h>
//#include <Rinternals.h> // included by siena07models.h
#include <R_ext/Rdynload.h>

extern "C"
{

/*
  From tools::package_native_routine_registration_skeleton.

  And I looked at how it is done in base package methods.
*/

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
   CALLDEF(Behavior, 2),
   CALLDEF(Continuous, 2),
   CALLDEF(Bipartite, 2),
   CALLDEF(ChangingCovariates, 2),
   CALLDEF(ChangingDyadicCovariates, 2),
   CALLDEF(clearStoredChains, 3),
   CALLDEF(ConstantCovariates, 2),
   CALLDEF(Constraints, 7),
   CALLDEF(deleteData, 1),
   CALLDEF(deleteModel, 1),
   CALLDEF(DyadicCovariates, 2),
   CALLDEF(effects, 2),
   CALLDEF(ExogEvent, 2),
   CALLDEF(forwardModel, 16),
   CALLDEF(getChainProbabilities, 8),
   CALLDEF(getTargets, 6),
   CALLDEF(interactionEffects, 2),
   CALLDEF(mlInitializeSubProcesses, 10),
   CALLDEF(mlMakeChains, 9),
   CALLDEF(mlPeriod, 14),
   CALLDEF(OneMode, 2),
   CALLDEF(setupData, 2),
   CALLDEF(setupModelOptions, 12),
    {NULL, NULL, 0}
};

void R_init_RSiena(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
