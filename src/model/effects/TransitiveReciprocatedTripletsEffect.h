/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveReciprocatedTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveReciprocatedTripletsEffect.
 *****************************************************************************/

#ifndef TRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_
#define TRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the reciprocated transitive triplets effect for non-symmetric
 * networks.
 */
class TransitiveReciprocatedTripletsEffect : public NetworkEffect
{
public:
	TransitiveReciprocatedTripletsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*TRANSITIVERECIPROCATEDTRIPLETSEFFECT_H_*/
