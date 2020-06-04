/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateMaximumBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * AltersCovariateMaximumEffect class.
 *****************************************************************************/

#ifndef ALTERSCOVARIATEMAXIMUMEFFECT_H_
#define ALTERSCOVARIATEMAXIMUMEFFECT_H_

#include "CovariateAndNetworkBehaviorEffect.h"

namespace siena
{


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Alters covariate maximum effect (see manual)
 *
 */
    class AltersCovariateMaximumEffect :
            public CovariateAndNetworkBehaviorEffect
    {
    public:
        AltersCovariateMaximumEffect(const EffectInfo * pEffectInfo);

        virtual double calculateChangeContribution(int actor,
                                                   int difference);
        virtual double egoEndowmentStatistic(int ego, const int * difference,
                                             double * currentValues);
        virtual double egoStatistic(int ego, double * currentValues);

    protected:

    private:

    };

}

#endif /*ALTERSCOVARIATEMAXIMUMEFFECT_H_*/
