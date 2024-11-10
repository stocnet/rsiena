/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: VarianceAlterEffect.h
 *
 * Description: This file contains the definition of the
 * VarianceAlterEffect class.
 *****************************************************************************/

#ifndef VARIANCEALTEREFFECT_H_
#define VARIANCEALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Variance alter effect defined as the product of the ego with the variance
 * of its neighbors (with respect to a certain network).
 */
class VarianceAlterEffect : public NetworkDependentBehaviorEffect
{
public:
    VarianceAlterEffect(const EffectInfo * pEffectInfo);

    virtual double calculateChangeContribution(int actor,
        int difference);
    virtual double egoEndowmentStatistic(int ego, const int * difference,
        double * currentValues);
    virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*VARIANCEALTEREFFECT_H_*/
