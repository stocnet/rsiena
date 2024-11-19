/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: VarianceAlterSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * VarianceAlterSimilarityEffect class.
 *****************************************************************************/
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "VarianceAlterSimilarityEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"


using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] average indicates if one of the average effects is required
 * @param[in] alterPopularity indicates if the similarity scores have to
 * be multiplied by the in-degrees of alters
 */
VarianceAlterSimilarityEffect::VarianceAlterSimilarityEffect(
    const EffectInfo * pEffectInfo,
    bool average,
    bool alterPopularity,
    bool egoPopularity) :
        NetworkDependentBehaviorEffect(pEffectInfo)
{
    this->laverage = average;
    this->lalterPopularity = alterPopularity;
    this->legoPopularity = egoPopularity;
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double VarianceAlterSimilarityEffect::calculateChangeContribution(int actor,
    int difference)
{
    double contribution = 0;
    double varianceAlt = 0;
    const Network * pNetwork = this->pNetwork();
    int neighborCount = pNetwork->outDegree(actor);
    
    if (neighborCount > 0)
    {
        // The formula for the interaction effect of
        // average similarity and variance alters is:
        // s_i(x) = [var(v_j) over all neighbors j of i] *
        // avg(sim(v_i, v_j) - centeringConstant) over all neighbors
        // j of i.
        // sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
        // We need to calculate the change delta in s_i(x), if we changed
        // v_i to v_i + d (d being the given amount of change in v_i).
        // To this end, we can disregard the centering constant and
        // compute the average change in similarity, namely,
        // avg(sim(v_i + d, v_j) - sim(v_i, v_j)) =
        // avg(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range) =
        // avg(|v_i-v_j| - |v_i+d-v_j|) / range,
        // the average being taken over all neighbors of i.
        // The reasoning for avg. similarity x variance alter effect is
        // similar.
        // This is what is calculated below.
        
        // step 1: determine the average similarity contribution
        int oldValue = this->value(actor);
        int newValue = oldValue + difference;
        int totalChange = 0;
        
        for (IncidentTieIterator iter = pNetwork->outTies(actor);
             iter.valid();
             iter.next())
        {
            int j = iter.actor();
            int alterValue = this->value(j);
            int change =
            std::abs(oldValue - alterValue) - std::abs(newValue - alterValue);
            
            if (this->lalterPopularity)
            {
                change *= pNetwork->inDegree(j);
            }
            
            totalChange += change;
        }
        
        contribution = ((double) totalChange) / this->range();
        
        if (this->laverage)
        {
            contribution /= pNetwork->outDegree(actor);
        }
        
        if (this->legoPopularity)
        {
            contribution *= pNetwork->inDegree(actor);
        }
        
        // step 2: compute the variance of the actor's alters
        for (IncidentTieIterator iter = pNetwork->outTies(actor);
             iter.valid();
             iter.next())
        {
            int j = iter.actor();
            varianceAlt += pow(this->centeredValue(j), 2);
        }
        varianceAlt /= neighborCount;
        varianceAlt -= pow(totalAlterValue(actor) / neighborCount, 2);
        varianceAlt -= this->variance();
    }
    
    // step 3: determine the effect contribution
    contribution *= varianceAlt;
    
    return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double VarianceAlterSimilarityEffect::egoStatistic(int ego,
    double * currentValues)
{
    const Network * pNetwork = this->pNetwork();

    double statistic = 0;
    double varianceAlt = 0;
    int neighborCount = pNetwork->outDegree(ego);

    if (neighborCount > 0)
    {
        double sum = 0, sumSquares = 0;

        for (IncidentTieIterator iter = pNetwork->outTies(ego);
             iter.valid();
             iter.next())
        {
            int j = iter.actor();
            sum += currentValues[j];
            sumSquares += pow(currentValues[j], 2);
            
            if (!this->missing(this->period(), j) &&
                !this->missing(this->period() + 1, j))
            {
                double tieStatistic =
                this->similarity(currentValues[ego], currentValues[j]);
                
                if (this->lalterPopularity)
                {
                    tieStatistic *= pNetwork->inDegree(j);
                }
                
                statistic += tieStatistic;
            }
        }
        
        varianceAlt = sumSquares / neighborCount - pow(sum / neighborCount, 2);
        
        if (this->laverage) // neighborCount > 0
        {
            statistic /= neighborCount;
        }
        
        if (this->legoPopularity)
        {
            statistic *= pNetwork->inDegree(ego);
        }
        
        statistic *= (varianceAlt - this->variance());
    }
    
    return statistic;
}


/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double VarianceAlterSimilarityEffect::egoEndowmentStatistic(int ego, const int * difference,
    double * currentValues)
{
    throw runtime_error(string("endowmentStatistic not implemented for") +
               "variance alter x average similarity.");
    // This egoEndowmentStatistic function can likely be removed.
}

}
