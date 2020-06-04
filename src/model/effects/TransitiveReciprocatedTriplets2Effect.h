/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveReciprocatedTriplets2Effect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveReciprocatedTriplets2Effect.
 *****************************************************************************/

#ifndef TRANSITIVERECIPROCATEDTRIPLETS2EFFECT_H_
#define TRANSITIVERECIPROCATEDTRIPLETS2EFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the reciprocated transitive triplets effect for non-symmetric
 * networks.
 */
class TransitiveReciprocatedTriplets2Effect : public NetworkEffect
{
public:
	TransitiveReciprocatedTriplets2Effect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*TRANSITIVERECIPROCATEDTRIPLETS2EFFECT_H_*/
