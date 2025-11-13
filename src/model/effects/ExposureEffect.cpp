#include <cmath>
#include <cstring>
#include "ExposureEffect.h"
#include "utils/Utils.h"
#include "model/variables/BehaviorVariable.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

namespace siena {

double ExposureEffect::proximityValue(const Network* pNetwork, int i, int period) const
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
            double alterValue = this->lpBehaviorVariable->value(iter.actor());
            if (alterValue >= 0.5)
            {
                numInfectedAlter++;
            }
            totalAlterValue += alterValue;
        }
    }

    totalAlterValue = this->applyInternalEffectParameter(totalAlterValue, numInfectedAlter);

    totalAlterValue *= egoNumer;

    double rawStatistic = (egoDenom > 1) ? totalAlterValue / egoDenom : totalAlterValue;
    return rawStatistic;
}

}