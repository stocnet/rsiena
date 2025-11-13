#include <cmath>
#include <cstring>
#include <stdexcept>
#include "InfectEffect.h"
#include "utils/Utils.h"
#include "model/variables/BehaviorVariable.h"
#include "network/OneModeNetwork.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

namespace siena {

double InfectEffect::proximityValue(const Network* pNetwork, int i, int period) const
{
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
            // infectIn: weighted by in-degree
            if (this->leffectName == "infectIn")
            {
                alterValue *= pNetwork->inDegree(iter.actor());
            }
            // infectDeg or infectOut: weighted by out-degree
            else if (this->leffectName == "infectDeg" || this->leffectName == "infectOut")
            {
                alterValue *= pNetwork->outDegree(iter.actor());
            }
             else if (this->leffectName == "infectCovar")
            {
                if (this->lpConstantCovariate)
                {
                    alterValue *= this->lpConstantCovariate->value(iter.actor());
                }
                else if (this->lpChangingCovariate)
                {
                    alterValue *= this->lpChangingCovariate->value(iter.actor(), period);
                }
                else
                {
                    throw std::logic_error("No individual covariate found for infectCovar.");
                }
            }
            totalAlterValue += alterValue;
        }
    }

    totalAlterValue = this->applyInternalEffectParameter(totalAlterValue, numInfectedAlter);

    return totalAlterValue;
}

}