#include "DiffusionRateEffectFactory.h"

namespace siena {

DiffusionRateEffect* DiffusionRateEffectFactory::create(
    const EffectInfo* pEffectInfo,
    const Network* pNetwork,
    const BehaviorVariable* pBehaviorVariable,
    const ConstantCovariate* pConstantCovariate,
    const ChangingCovariate* pChangingCovariate,
    const std::string& effectName,
    double parameter,
    double internalEffectParameter)
{
    if (effectName == "avExposure" || effectName == "totExposure")
    {
        return new ExposureEffect(pEffectInfo, pNetwork, pBehaviorVariable, effectName, parameter, internalEffectParameter);
    }
    else if (effectName == "infectIn" || effectName == "infectDeg" || effectName == "infectOut" || effectName == "infectCovar")
    {
        return new InfectEffect(pEffectInfo, pNetwork, pBehaviorVariable, effectName, parameter, internalEffectParameter);
    }
    else if (effectName == "susceptAvIn" || effectName == "susceptAvCovar")
    {
        return new SusceptibilityEffect(pEffectInfo, pNetwork, pBehaviorVariable, pConstantCovariate, pChangingCovariate, effectName, parameter, internalEffectParameter);
    }
    else if (effectName == "anyInExposureDist2" || effectName == "totInExposureDist2" ||
             effectName == "avTinExposureDist2" || effectName == "totAInExposureDist2")
    {
        return new Distance2ExposureEffect(pEffectInfo, pNetwork, pBehaviorVariable, effectName, parameter, internalEffectParameter);
    }
    else
    {
        throw std::domain_error("Unknown diffusion rate effect: " + effectName);
    }
}

}