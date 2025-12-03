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
#include <R_ext/Print.h>

namespace siena {

double ExposureEffect::proximityValue(const Network* pNetwork, int i) const
{
    int egoNumer = 1;
    int egoDenom = 1;

    if (this->leffectName == "avExposure")
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

    totalAlterValue = this->applyThreshold(totalAlterValue, numInfectedAlter);

    totalAlterValue *= egoNumer;

    double rawStatistic = (egoDenom > 1) ? totalAlterValue / egoDenom : totalAlterValue;
    return rawStatistic;
}

}