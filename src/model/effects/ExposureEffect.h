#ifndef EXPOSUREEFFECT_H_
#define EXPOSUREEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class ExposureEffect : public DiffusionRateEffect
{
public:
    using DiffusionRateEffect::DiffusionRateEffect;
protected:
    double proximityValue(const Network* pNetwork, int i) const;
};

}

#endif /* EXPOSUREEFFECT_H_ */