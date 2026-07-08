/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ExposureEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * ExposureEffect, which handles direct exposure-based diffusion rate effects.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include "ExposureEffect.h"
#include "utils/Utils.h"
#include "model/variables/BehaviorVariable.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include <R_ext/Print.h>

namespace siena {
ExposureEffect::ExposureEffect(const EffectInfo* pEffectInfo)
    : DiffusionRateEffect(pEffectInfo)
{
    // 1. Process threshold configuration once
    int rawParam = static_cast<int>(std::round(pEffectInfo->internalEffectParameter()));
    this->labsThreshold = std::abs(rawParam);
    this->lhasThreshold = (rawParam != 0);
    this->lcapAtThreshold = (rawParam < 0);

    // 2. Process string match once. The heavy character checking ends here!
    this->lIsAvExposure = (this->leffectName == "avExposure");
}

/**
 * Wrapper to fix abstract threshold utility to concrete threshold values.
 */
double ExposureEffect::applyThreshold(double value, int numInfectedAlter) const
{
    return DiffusionRateEffect::applyThreshold(value, numInfectedAlter, this->labsThreshold, this->lcapAtThreshold);
}

double ExposureEffect::proximityValue(const Network* pNetwork, int i) const
{
    int egoNumer = 1;
    int egoDenom = 1;

    if (this->lIsAvExposure)
    {
        egoDenom = std::max(1, pNetwork->outDegree(i));
    }

    double totalAlterValue = 0;
    int numInfectedAlter = 0;

    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
        {
            double alterValue = this->value(iter.actor());
            if (alterValue >= 0.5)
            {
                numInfectedAlter++;
            }
            totalAlterValue += alterValue;
        }
    }

    if (this->lhasThreshold)
    {
        totalAlterValue = this->applyThreshold(totalAlterValue, numInfectedAlter);
    }
    totalAlterValue *= egoNumer;

    return (egoDenom > 1) ? totalAlterValue / egoDenom : 
                            totalAlterValue;
}

}