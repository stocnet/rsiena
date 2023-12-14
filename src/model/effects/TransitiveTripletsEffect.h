/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveTripletsEffect.
 *****************************************************************************/

#ifndef TRANSITIVETRIPLETSEFFECT_H_
#define TRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive triplets effect for non-symmetric
 * networks.
 */
class TransitiveTripletsEffect : public NetworkEffect
{
public:
	TransitiveTripletsEffect(const EffectInfo * pEffectInfo, 
					bool twoPath, bool twoInStar);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	
private:
	// twoPath indicates whether twopaths will be closed,
	// twoInStar indicates whether two-instars will be closed.
	bool ltwoPath {};
	bool ltwoInStar {};
};

}

#endif /*TRANSITIVETRIPLETSEFFECT_H_*/
