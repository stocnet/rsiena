/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * SimilarityEffect class.
 *****************************************************************************/

#ifndef VARIANCEALTERSIMILARITYEFFECT_H_
#define VARIANCEALTERSIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * This class could implement several behavior effects related to the interaction
 * between similarity and the variance alter effect. Currently, it only implements
 * the average simiarity x variance alter effect. However, similar to the regular
 * avSim effect, this could also implement (see manual):
 * - Average similarity x popularity alter x variance alter
 * - Total similarity x variance alter
 * - Total similarity x popularity alter x variance alter
 * - Average similarity x popularity ego x variance alter
 */
class VarianceAlterSimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
    VarianceAlterSimilarityEffect(const EffectInfo * pEffectInfo,
        bool average,
        bool alterPopularity,
        bool egoPopularity);

    virtual double calculateChangeContribution(int actor,
        int difference);
    virtual double egoEndowmentStatistic(int ego, const int * difference,
        double * currentValues);
    virtual double egoStatistic(int ego, double * currentValues);

private:
    bool laverage;
    bool lalterPopularity;
    bool legoPopularity;
};

}

#endif /*VARIANCEALTERSIMILARITYEFFECT_H_*/
