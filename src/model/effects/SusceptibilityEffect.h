/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SusceptibilityEffect.h
 *
 * Description: This file contains the definition of the class
 * SusceptibilityEffect, which handles susceptibility diffusion rate effects.
 *****************************************************************************/

#ifndef SUSCEPTIBILITYEFFECT_H_
#define SUSCEPTIBILITYEFFECT_H_

#include "DiffusionRateEffect.h"

namespace siena
{

/**
 * Handles all susceptibility-based diffusion rate effects (e.g., susceptAvIn, susceptAvCovar).
 */
class SusceptibilityEffect : public DiffusionRateEffect
{
public:
    using DiffusionRateEffect::DiffusionRateEffect; // Inherit constructors etc.

protected:
    double proximityValue(const Network* pNetwork, int i) const;
};

}

#endif /* SUSCEPTIBILITYEFFECT_H_ */