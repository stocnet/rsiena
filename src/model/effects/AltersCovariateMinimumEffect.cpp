/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateMinimumEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateMinimumEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "AltersCovariateMinimumEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
    AltersCovariateMinimumEffect::AltersCovariateMinimumEffect(
            const EffectInfo * pEffectInfo) :
            CovariateAndNetworkBehaviorEffect(pEffectInfo)
    {

    }



/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
    double AltersCovariateMinimumEffect::calculateChangeContribution(int actor,
                                                                     int difference)
    {

        double statistic = difference * this->minimumAlterValue(actor);

        return statistic;
    }

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
    double AltersCovariateMinimumEffect::egoStatistic(int ego, double * currentValues)
    {
        double statistic = currentValues[ego] * this->minimumAlterValue(ego);

        return statistic;
    }

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
    double AltersCovariateMinimumEffect::egoEndowmentStatistic(int ego,
                                                               const int * difference,
                                                               double * currentValues)
    {
        double statistic = 0;

        if (difference[ego] > 0 && !this->missingDummy(ego))
        {
                statistic -= difference[ego] * this->minimumAlterValue(ego);

        }
        return statistic;
    }


}
