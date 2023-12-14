/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InInDegreeAssortativityEffect.h
 *
 * Description: This file contains the definition of the
 * InInDegreeAssortativityEffect class.
 *****************************************************************************/

#ifndef ININDEGREEASSORTATIVITYEFFECT_H_
#define ININDEGREEASSORTATIVITYEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The in-in degree assortativity effect (see manual).
 */
class InInDegreeAssortativityEffect : public NetworkEffect
{
public:
	InInDegreeAssortativityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of degrees must be used
	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*ININDEGREEASSORTATIVITYEFFECT_H_*/
