/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DensityEffect.h
 *
 * Description: This file contains the declaration of the class DensityEffect.
 *****************************************************************************/

#ifndef DENSITYEFFECT_H_
#define DENSITYEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the density effect.
 */
class DensityEffect : public NetworkEffect
{
public:
	DensityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
	virtual bool egoEffect() const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*DENSITYEFFECT_H_*/
