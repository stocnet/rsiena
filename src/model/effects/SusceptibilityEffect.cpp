#include <cmath>
#include <cstring>
#include <stdexcept>
#include "SusceptibilityEffect.h"
#include "utils/Utils.h"
#include "model/variables/BehaviorVariable.h"
#include "network/OneModeNetwork.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

namespace siena
{

double SusceptibilityEffect::proximityValue(const Network* pNetwork, int i, int period) const
{
    int egoNumer = 1;
    int egoDenom = 1;

    // susceptAvIn: numerator is in-degree, denominator is out-degree
    if (this->leffectName == "susceptAvIn")
    {
        egoNumer = pNetwork->inDegree(i);
        egoDenom = std::max(1, pNetwork->outDegree(i));
    }
    // susceptAvCovar: denominator is out-degree
    else if (this->leffectName == "susceptAvCovar")
    {
        egoDenom = std::max(1, pNetwork->outDegree(i));
    }

    double totalAlterValue = 0;
    int numInfectedAlter = 0;

    if (pNetwork->outDegree(i) > 0)
    {
        for (IncidentTieIterator iter = pNetwork->outTies(i); iter.valid(); iter.next())
        {
            double alterValue = this->lpBehaviorVariable->value(iter.actor());
            if (alterValue >= 0.5)
            {
                numInfectedAlter++;
            }
            totalAlterValue += alterValue;
        }
    }

    // Internal effect parameter thresholding
    totalAlterValue = this->applyInternalEffectParameter(totalAlterValue, numInfectedAlter);

    totalAlterValue *= egoNumer;

    double rawStatistic;
    if (egoDenom > 1)
    {
        rawStatistic = totalAlterValue / egoDenom;
    }
    else
    {
        rawStatistic = totalAlterValue;
    }

    if (this->leffectName == "susceptAvCovar")
    {
        if (this->lpConstantCovariate)
        {
            rawStatistic *= this->lpConstantCovariate->value(i);
        }
        else if (this->lpChangingCovariate)
        {
            rawStatistic *= this->lpChangingCovariate->value(i, period);
        }
        else
        {
            throw std::logic_error("No individual covariate found for susceptAvCovar.");
        }
    }

    return rawStatistic;
}

}