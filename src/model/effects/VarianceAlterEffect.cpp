/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: VarianceAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * VarianceAlterEffect class.
 *****************************************************************************/

#include <R_ext/Print.h>
#include <cmath>
#include "VarianceAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{
    /**
     * Constructor.
     */
    VarianceAlterEffect::VarianceAlterEffect(const EffectInfo * pEffectInfo) :
    NetworkDependentBehaviorEffect(pEffectInfo)
    {

    }
    
    
    /**
     * Calculates the change in the statistic corresponding to this effect if
     * the given actor would change his behavior by the given amount.
     */
    double VarianceAlterEffect::calculateChangeContribution(int actor,
                                                           int difference)
    {
        double contribution = 0;
        const Network * pNetwork = this->pNetwork();
        int neighborCount = pNetwork->outDegree(actor);

        if (neighborCount > 1) // for 1 neighbor, variance is 0
        {
            // The formula for the effect:
            // s_i(x) = v_i * ( [var(v_j) over all neighbors j of i] - overall var ).
            // The overall variance serves as centering constant.
            // We need to calculate the change delta in s_i(x), if we changed
            // v_i to v_i + d (d being the given amount of change in v_i).
            // This is d * (var(v_j) - overall var),
            // the first variance being taken over all neighbors of i.
            // This is what is calculated below.
            for (IncidentTieIterator iter = pNetwork->outTies(actor);
                     iter.valid();
                     iter.next())
            {
                int j = iter.actor();
                contribution += pow(this->centeredValue(j), 2);
            }
            contribution /= neighborCount;
            contribution -= pow(totalAlterValue(actor) / neighborCount, 2);
        }

        if (neighborCount > 0)
        {
            contribution -= this->variance();
            contribution *= difference;
        }
        return contribution;
    }
    
    
    /**
     * Returns the statistic corresponding to the given ego with respect to the
     * given values of the behavior variable.
     */
    double VarianceAlterEffect::egoStatistic(int i, double * currentValues)
    {
        const Network * pNetwork = this->pNetwork();
        int neighborCount = pNetwork->outDegree(i);
        double statistic = 0;
     //   Rprintf("--------------- Actor %d ----------------\n", i);
      //  Rprintf("Current value: %f\n", currentValues[i]);
        if (neighborCount > 1)
        {
            double sum = 0, sumSquares = 0;
        
            for (IncidentTieIterator iter = pNetwork->outTies(i);
                iter.valid();
                iter.next())
            {
                int j = iter.actor();
                sum += currentValues[j];
                sumSquares += pow(currentValues[j], 2);
        //        Rprintf("%d is a neighbor with score %f\n", j, currentValues[j]);
            }

            statistic = sumSquares / neighborCount - pow(sum / neighborCount, 2);
        }
        
        if (neighborCount > 0)
        {
            statistic -= this->variance();
            statistic *= currentValues[i];
        }
        
        return statistic;
    }
    
    /**
     * Returns the statistic corresponding to the given ego as part of
     * the endowment function with respect to the initial values of a
     * behavior variable and the current values.
     */
    double VarianceAlterEffect::egoEndowmentStatistic(int ego,
                                                     const int * difference,
                                                     double * currentValues)
    {
        double statistic = 0;
        const Network * pNetwork = this->pNetwork();
        int neighborCount = pNetwork->outDegree(ego);

        if (difference[ego] > 0 && neighborCount > 1)
        {
            double thisStatistic = 0;
            double previousStatistic = 0;
            double thisSum = 0, thisSumSquares = 0;
            double previousSum = 0, previousSumSquares = 0;
            
            for (IncidentTieIterator iter = pNetwork->outTies(ego);
                 iter.valid();
                 iter.next())
            {
                int j = iter.actor();
                
                double alterValue = currentValues[j];
                double alterPreviousValue = currentValues[j] + difference[j];
                
                thisSum += alterValue;
                thisSumSquares += pow(alterValue, 2);
                previousSum += alterPreviousValue;
                previousSumSquares += pow(alterPreviousValue, 2);
            }
            
            thisStatistic = thisSumSquares / neighborCount -
                pow(thisSum / neighborCount, 2);
            thisStatistic *= currentValues[ego];
            
            previousStatistic = previousSumSquares / neighborCount -
                pow(previousSum / neighborCount, 2);
            previousStatistic *= (currentValues[ego] + difference[ego]);
            
            statistic = thisStatistic - previousStatistic;
        }
        
        if (neighborCount > 0)
        {
            statistic += (difference[ego] * this->variance());
        }
        return statistic;
    }

}
