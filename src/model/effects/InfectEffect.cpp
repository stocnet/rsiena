/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InfectEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * InfectEffect, which handles infection-based diffusion rate effects.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include <stdexcept>
#include "InfectEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include <R_ext/Print.h>

namespace siena {

InfectEffect::InfectEffect(const EffectInfo* pEffectInfo)
    : DiffusionRateEffect(pEffectInfo)
{
    int rawParam = static_cast<int>(std::round(pEffectInfo->internalEffectParameter()));
    this->labsThreshold = std::abs(rawParam);
    this->lhasThreshold = (rawParam != 0);
    this->lcapAtThreshold = (rawParam < 0);

    // Cache string identifiers as lightweight booleans
    this->lIsInfectIn = (this->leffectName == "infectIn");
    this->lIsInfectDegOrOut = (this->leffectName == "infectDeg" || this->leffectName == "infectOut");
    this->lIsInfectCovar = (this->leffectName == "infectCovar");
}

void InfectEffect::initialize(const Data* pData, State* pState, int period, Cache* pCache)
{
    // 1. Parent binds the network pointers and covariates
    DiffusionRateEffect::initialize(pData, pState, period, pCache);

    // 2. Safely validate covariate presence once per wave/period setup
    if (this->lIsInfectCovar && !this->lpConstantCovariate && !this->lpChangingCovariate)
    {
        throw std::logic_error("Covariate data missing for infectCovar.");
    }
}

double InfectEffect::applyThreshold(double value, int numInfectedAlter) const
{
    return DiffusionRateEffect::applyThreshold(value, numInfectedAlter, this->labsThreshold, this->lcapAtThreshold);
}

double InfectEffect::proximityValue(const Network* pNetwork, int i)
{
    double totalAlterValue = 0;
    int numInfectedAlter = 0;

    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
        {
            int alterActor = iter.actor();
            double alterValue = this->value(alterActor);

            if (alterValue >= 0.5)
            {
                numInfectedAlter++;
            }

            if (this->lIsInfectIn)
            {
                alterValue *= pNetwork->inDegree(alterActor);
            }
            else if (this->lIsInfectDegOrOut)
            {
                alterValue *= pNetwork->outDegree(alterActor);
            }
            else if (this->lIsInfectCovar)
            {
                alterValue *= this->lpChangingCovariate ? 
                              this->lpChangingCovariate->value(alterActor, this->period()) : 
                              this->lpConstantCovariate->value(alterActor);
            }

            totalAlterValue += alterValue;
        }
    }

    if (this->lhasThreshold)
    {
        totalAlterValue = this->applyThreshold(totalAlterValue, numInfectedAlter);
    }

    return totalAlterValue;
}

}