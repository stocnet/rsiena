/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Distance2ExposureEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * Distance2ExposureEffect, which handles distance-2 exposure rate effects.
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include "Distance2ExposureEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include <R_ext/Print.h>

namespace siena
{

Distance2ExposureEffect::Distance2ExposureEffect(const EffectInfo* pEffectInfo) :
    DiffusionRateEffect(pEffectInfo)
{
    // 1. Process threshold configuration once
    int rawParam = static_cast<int>(std::round(pEffectInfo->internalEffectParameter()));
    this->labsThreshold = std::abs(rawParam);
    this->lhasThreshold = (rawParam != 0);
    this->lcapAtThreshold = (rawParam < 0);

    // 2. Clear out all hot-path string comparisons by mapping them to flags
    this->lIsAvTinExposureDist2 = (this->leffectName == "avTinExposureDist2");
    this->lIsTotAInExposureDist2 = (this->leffectName == "totAInExposureDist2");
    this->lIsAnyInExposureDist2 = (this->leffectName == "anyInExposureDist2");
}

void Distance2ExposureEffect::initialize(const Data* pData, State* pState, int period, Cache* pCache)
{
    DiffusionRateEffect::initialize(pData, pState, period, pCache);
}

double Distance2ExposureEffect::applyThreshold(double value, int numInfectedAlter) const
{
    return DiffusionRateEffect::applyThreshold(value, numInfectedAlter, this->labsThreshold, this->lcapAtThreshold);
}

double Distance2ExposureEffect::proximityValue(const Network* pNetwork, int i) const
{
    int egoDenom = 1;

    if (this->lIsAvTinExposureDist2)
    {
        egoDenom = std::max(1, pNetwork->outDegree(i));
    }

    double totalAlterValue = 0;
    int numInfectedAlter = 0;

    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
        {
            int j = iter.actor();
            double totalAlterInDist2Value = 0;

            for (IncidentTieIterator iterH = pNetwork->inTies(j); iterH.valid(); iterH.next())
            {
                int alterH = iterH.actor();
                if (i != alterH)
                {
                    double alterInDist2Value = this->value(alterH);
                    if (alterInDist2Value >= 0.5)
                    {
                        numInfectedAlter++;
                    }
                    totalAlterInDist2Value += alterInDist2Value;
                }
            }

            if (this->lIsTotAInExposureDist2)
            {
                int normalizedInDegree = pNetwork->inDegree(j) - 1;
                if (normalizedInDegree > 0)
                {
                    totalAlterInDist2Value /= normalizedInDegree;
                }
            }
            else if (this->lIsAnyInExposureDist2)
            {
                totalAlterInDist2Value = std::min(totalAlterInDist2Value, 1.0);
            }

            totalAlterValue += totalAlterInDist2Value;
        }
    }

    if (this->lhasThreshold)
    {
        totalAlterValue = this->applyThreshold(totalAlterValue, numInfectedAlter);
    }

    return (egoDenom > 1) ? totalAlterValue / egoDenom : totalAlterValue;
}

}