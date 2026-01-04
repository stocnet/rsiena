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
#include <R_ext/Print.h>

namespace siena
{

double Distance2ExposureEffect::proximityValue(const Network* pNetwork, int i) const
{
    int egoDenom = 1;

    // avTinExposureDist2 is normalized by outdegree
    if (this->leffectName == "avTinExposureDist2")
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
                if (i != iterH.actor())
                {
                    double alterInDist2Value = this->value(iterH.actor());
                    if (alterInDist2Value >= 0.5)
                    {
                        numInfectedAlter++;
                    }
                    totalAlterInDist2Value += alterInDist2Value;
                }
            }

            // totAInExposureDist2: average over in-degree minus self
            if ((this->leffectName == "totAInExposureDist2") && ((pNetwork->inDegree(j) - 1) > 0))
            {
                totalAlterInDist2Value /= (pNetwork->inDegree(j) - 1);
            }
            // anyInExposureDist2: at least one infected in-neighbor
            if (this->leffectName == "anyInExposureDist2")
            {
                totalAlterInDist2Value = std::min(totalAlterInDist2Value, 1.0);
            }

            totalAlterValue += totalAlterInDist2Value;
        }
    }

    // does using the thresholds make sense here?
    totalAlterValue = this->applyThreshold(totalAlterValue, numInfectedAlter);

    double rawStatistic;
    if (egoDenom > 1)
    {
        rawStatistic = totalAlterValue / egoDenom;
    }
    else
    {
        rawStatistic = totalAlterValue;
    }

    return rawStatistic;
}

}