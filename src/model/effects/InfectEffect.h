/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InfectEffect.h
 *
 * Description: This file contains the definition of the class
 * InfectEffect, which handles infection-based diffusion rate effects.
 *****************************************************************************/

#ifndef INFECTEFFECT_H_
#define INFECTEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena {

class InfectEffect : public DiffusionRateEffect
{
public:
    using DiffusionRateEffect::DiffusionRateEffect;

protected:
    double proximityValue(const Network* pNetwork, int i) const;
};

}

#endif /* INFECTEFFECT_H_ */