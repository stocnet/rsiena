#ifndef DISTANCE2EXPOSUREEFFECT_H_
#define DISTANCE2EXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena
{

/**
 * Handles all distance-2 diffusion rate effects (e.g., anyInExposureDist2, totInExposureDist2, avTinExposureDist2, totAInExposureDist2).
 */
class Distance2ExposureEffect : public DiffusionRateEffect
{
public:
    using DiffusionRateEffect::DiffusionRateEffect; // Inherit constructors

protected:
    double proximityValue(const Network* pNetwork, int i, int period) const;
};

}

#endif /* DISTANCE2EXPOSUREEFFECT_H_ */