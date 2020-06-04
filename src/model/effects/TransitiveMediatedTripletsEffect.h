/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveMediatedTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveMediatedTripletsEffect.
 *****************************************************************************/

#ifndef TRANSITIVEMEDIATEDTRIPLETSEFFECT_H_
#define TRANSITIVEMEDIATEDTRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive mediated triplets effect.
 */
class TransitiveMediatedTripletsEffect : public NetworkEffect
{
public:
	TransitiveMediatedTripletsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*TRANSITIVEMEDIATEDTRIPLETSEFFECT_H_*/
