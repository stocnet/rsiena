/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ContinuousInteractionEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * ContinuousInteractionEffect.
 *****************************************************************************/

#include "ContinuousInteractionEffect.h"
#include <cmath>

namespace siena
{

/**
 * Constructs a new interaction effect between the given effects.
 * The parameter pEffect3 should be 0 for two-way interactions.
 * This effect takes the ownership of the given effects, which mean
 * that the given effects will be destroyed as soon as this
 * effect is destroyed.
 */
ContinuousInteractionEffect::ContinuousInteractionEffect(
                                                         const EffectInfo * pEffectInfo,
                                                         ContinuousEffect * pEffect1,
                                                         ContinuousEffect * pEffect2,
                                                         ContinuousEffect * pEffect3) : ContinuousEffect(pEffectInfo)
{
    this->lpEffect1 = pEffect1;
    this->lpEffect2 = pEffect2;
    this->lpEffect3 = pEffect3;
}


/**
 * Deallocates this effects.
 */
ContinuousInteractionEffect::~ContinuousInteractionEffect()
{
    delete this->lpEffect1;
    delete this->lpEffect2;
    delete this->lpEffect3;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ContinuousInteractionEffect::initialize(const Data * pData,
                                             State * pState,
                                             int period,
                                             Cache * pCache)
{
    ContinuousEffect::initialize(pData, pState, period, pCache);
    this->lpEffect1->initialize(pData, pState, period, pCache);
    this->lpEffect2->initialize(pData, pState, period, pCache);
    
    if (this->lpEffect3)
    {
        this->lpEffect3->initialize(pData, pState, period, pCache);
    }
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling BehaviorEffect::calculateChangeContribution(...).
 */
void ContinuousInteractionEffect::preprocessEgo(int ego)
{
    ContinuousEffect::preprocessEgo(ego);
    
    this->lpEffect1->preprocessEgo(ego);
    this->lpEffect2->preprocessEgo(ego);
    
    if (this->lpEffect3)
    {
        this->lpEffect3->preprocessEgo(ego);
    }
}


/**
 * Calculates the change in the statistic corresponding to this effecf if
 * the given actor were to change his behavior by the given amount.
 */
double ContinuousInteractionEffect::calculateChangeContribution(int actor)
{
    double contribution =
    this->lpEffect1->calculateChangeContribution(actor) *
    this->lpEffect2->calculateChangeContribution(actor);
    
    if (this->lpEffect3)
    {
        contribution *= this->lpEffect3->calculateChangeContribution(actor);
    }
    
    return contribution;
}


/**
 * Calculates the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double ContinuousInteractionEffect::egoStatistic(int i,
                                                 double * currentValues)
{
    double statistic;
    
    if (std::abs(currentValues[i]) > 1e-10)
    {
        statistic = this->lpEffect1->egoStatistic(i, currentValues) *
        this->lpEffect2->egoStatistic(i, currentValues);
        
        statistic /= currentValues[i];
        
        if (this->lpEffect3)
        {
            statistic *= this->lpEffect3->egoStatistic(i, currentValues)
            / currentValues[i];
        }
    }
    else
    {
        statistic = 0;
    }
    
    return statistic;
}

}
