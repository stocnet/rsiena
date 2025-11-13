#ifndef INFECTEFFECT_H_
#define INFECTEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class InfectEffect : public DiffusionRateEffect
{
public:
    using DiffusionRateEffect::DiffusionRateEffect;

protected:
    double proximityValue(const Network* pNetwork, int i, int period) const;
};

}

#endif /* INFECTEFFECT_H_ */