/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTriadsEffect.h
 *
 * Description: This file contains the declaration of the class
 * TransitiveTriadsEffect.
 *****************************************************************************/

#ifndef TRANSITIVETRIADSEFFECT_H_
#define TRANSITIVETRIADSEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class TransitiveTriadsEffect : public NetworkEffect
{
public:
	TransitiveTriadsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*TRANSITIVETRIADSEFFECT_H_*/
