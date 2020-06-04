/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateMinimumBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateMinimumEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEMINIMUMEFFECT_H_
#define ALTERSCOVARIATEMINIMUMEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate minimum effect (see manual)
 *
 */
    class AltersCovariateMinimumEffect :
            public CovariateAndNetworkBehaviorEffect
    {
    public:
        AltersCovariateMinimumEffect(const EffectInfo * pEffectInfo);

        virtual double calculateChangeContribution(int actor,
                                                   int difference);
        virtual double egoEndowmentStatistic(int ego, const int * difference,
                                             double * currentValues);
        virtual double egoStatistic(int ego, double * currentValues);

    protected:

    private:

    };

}

#endif /*ALTERSCOVARIATEMINIMUMEFFECT_H_*/
