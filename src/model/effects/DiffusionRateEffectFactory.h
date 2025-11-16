#ifndef DIFFUSIONRATEEFFECTFACTORY_H_
#define DIFFUSIONRATEEFFECTFACTORY_H_

#include "DiffusionRateEffect.h"
#include "ExposureEffect.h"
#include "InfectEffect.h"
#include "SusceptibilityEffect.h"
#include "Distance2ExposureEffect.h"

namespace siena {

class DiffusionRateEffectFactory
{
public:
    static DiffusionRateEffect* create(
        const EffectInfo* pEffectInfo,
        const Network* pNetwork,
        const BehaviorVariable* pBehaviorVariable,
        const ConstantCovariate* pConstantCovariate,
        const ChangingCovariate* pChangingCovariate,
        const std::string& effectName,
        double parameter,
        double internalEffectParameter);
};

}

#endif