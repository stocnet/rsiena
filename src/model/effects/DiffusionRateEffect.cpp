/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DiffusionRateEffect.
 *****************************************************************************/

#include <stdexcept>

#include <cmath>
#include <cstring>
#include "DiffusionRateEffect.h"
#include "utils/Utils.h"
#include "data/Data.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/tables/Cache.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
DiffusionRateEffect::DiffusionRateEffect(const EffectInfo* pEffectInfo)
    : BehaviorRateEffect(pEffectInfo)
{
    // Nothing else needed here; everything is set up in initialize
}

/**
 * Deallocates this effect object;
 */
DiffusionRateEffect::~DiffusionRateEffect()
{
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DiffusionRateEffect::initialize(
    const Data* pData, 
    State* pState, 
    int period, 
    Cache* pCache)
{
    BehaviorRateEffect::initialize(pData, pState, period, pCache);
    string networkName = this->pEffectInfo()->interactionName1();
    this->lpNetwork = pState->pNetwork(networkName);

    if (!this->lpNetwork) {
		throw logic_error("Network '" + networkName + "' expected.");
	}

    string name = this->pEffectInfo()->interactionName2();
    // could just be done for those effects that need it?
    if(!name.empty())
    {
        this->lpConstantCovariate = pData->pConstantCovariate(name);
	    this->lpChangingCovariate = pData->pChangingCovariate(name);
    }
    this->leffectName = this->pEffectInfo()->effectName(); // not always necessary?
    this->linternalEffectParameter = this->pEffectInfo()->internalEffectParameter();
    this->linternalEffectParameter = round(this->linternalEffectParameter); // why round?
    this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
    this->linternalNonZero = (this->linternalEffectParameter != 0);


}


// /**
//  * Constructor.
//  * @param[in] pNetwork the network this effect depends on
//  * @param[in] pBehaviorVariable the behavior variable this effect depends on
//  * @param[in] pConstantCovariate the covariate this effect depends on
//  * @param[in] pChangingCovariate the changing covariate this effect depends on
//  * @param[in] effectName the name of this effect
//  * @param[in] parameter the statistical parameter of this effect
//  * @param[in] internalEffectParameter the internal effect parameter
//  */
// DiffusionRateEffect::DiffusionRateEffect(const EffectInfo* pEffectInfo,
//     const Network * pNetwork,
//     const BehaviorVariable * pBehaviorVariable,
//     const ConstantCovariate * pConstantCovariate,
//     const ChangingCovariate * pChangingCovariate,
//     string effectName,
//     double parameter,
//     double internalEffectParameter) : Effect(pEffectInfo)
// {
//     this->lpNetwork = pNetwork;
//     this->lpBehaviorVariable = pBehaviorVariable;
//     this->lpChangingCovariate = pChangingCovariate;
//     this->lpConstantCovariate = pConstantCovariate;
//     this->leffectName = effectName;
//     this->lparameter = parameter;
//     this->linternalEffectParameter = round(internalEffectParameter);
//     this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
//     this->linternalNonZero = (this->linternalEffectParameter != 0);

//     if ((effectName == "infectCovar") && (this->linternalEffectParameter < 0))
//     {
//         throw logic_error("Negative internal parameter not permitted for effect "+effectName);
//     }
// }



double DiffusionRateEffect::proximityValue(const Network* pNetwork, int i) const
{
    throw std::logic_error("proximityValue not implemented for this effect type.");
}

double DiffusionRateEffect::applyThreshold(double value, int numInfectedAlter) const
{
    if (this->linternalNonZero)
    {
        if (numInfectedAlter < this->labsInternalEffectParameter)
        {
            value = 0;
        }
        else if (this->linternalEffectParameter < 0)
        {
            if (value > this->labsInternalEffectParameter)
            {
                value = this->labsInternalEffectParameter;
            }
        }
    }
    return value;
}

/**
 * Returns the rate contribution (for scores and statistics calculation).
 * Simple transformer: just returns the proximity value.
 */
double DiffusionRateEffect::calculateContribution(int i) const
{
    return this->proximityValue(this->lpNetwork, i);
}

/**
 * Stores the internal effect parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::setInternalEffectParameter(int parValue)
{
    this->linternalEffectParameter = parValue;
    this->labsInternalEffectParameter = std::abs(this->linternalEffectParameter);
    this->linternalNonZero = (this->linternalEffectParameter != 0);
}

/**
 * Returns the internal effect parameter for the diffusion rate effect.
 */
int DiffusionRateEffect::getInternalEffectParameter() const
{
    return this->linternalEffectParameter;
}

}