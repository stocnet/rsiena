/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SusceptibilityEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * SusceptibilityEffect, which handles susceptibility diffusion rate effects.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include <stdexcept>
#include "SusceptibilityEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include <R_ext/Print.h>

namespace siena
{

// --- CONSTRUCTOR: Parsed once at runtime birth ---
SusceptibilityEffect::SusceptibilityEffect(const EffectInfo* pEffectInfo)
    : DiffusionRateEffect(pEffectInfo)
{
    int rawParam = static_cast<int>(std::round(pEffectInfo->internalEffectParameter()));
    this->labsThreshold = std::abs(rawParam);
    this->lhasThreshold = (rawParam != 0);
    this->lcapAtThreshold = (rawParam < 0);

    this->lIsSusceptAvIn = (this->leffectName == "susceptAvIn");
    this->lIsSusceptAvCovar = (this->leffectName == "susceptAvCovar");

}

void SusceptibilityEffect::initialize(const Data* pData, State* pState, int period, Cache* pCache)
{
    DiffusionRateEffect::initialize(pData, pState, period, pCache);

    if (this->lIsSusceptAvCovar && !this->lpConstantCovariate && !this->lpChangingCovariate)
    {
        throw std::logic_error("Covariate data missing for susceptAvCovar.");
    }
}

// --- HOT PATH WRAPPER ---
double SusceptibilityEffect::applyThreshold(double value, int numInfectedAlter) const
{
    return DiffusionRateEffect::applyThreshold(value, numInfectedAlter, this->labsThreshold, this->lcapAtThreshold);
}

// --- SIMULATION EXECUTION: Executed millions of times ---
double SusceptibilityEffect::proximityValue(const Network* pNetwork, int i) const
{
    int egoNumer = 1;
    int egoDenom = 1;

    if (this->lIsSusceptAvIn)
    {
        egoNumer = pNetwork->inDegree(i);
        egoDenom = std::max(1, pNetwork->outDegree(i));
    }
    else if (this->lIsSusceptAvCovar)
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

    double rawStatistic = (egoDenom > 1) ? totalAlterValue / egoDenom : totalAlterValue;

    if (this->lIsSusceptAvCovar)
        {
            // Safe to execute without checking for total nulls because the constructor guaranteed validity!
            rawStatistic *= this->lpChangingCovariate ? 
                            this->lpChangingCovariate->value(i, this->period()) : 
                            this->lpConstantCovariate->value(i);
        }

    return rawStatistic;
}

}