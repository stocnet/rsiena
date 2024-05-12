/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InOutDegreeAssortativityEffect.h
 *
 * Description: This file contains the definition of the
 * InOutDegreeAssortativityEffect class.
 *****************************************************************************/

#ifndef INOUTDEGREEASSORTATIVITYEFFECT_H_
#define INOUTDEGREEASSORTATIVITYEFFECT_H_

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
 * The in-out degree assortativity effect (see manual).
 */
class InOutDegreeAssortativityEffect : public NetworkEffect
{
public:
	InOutDegreeAssortativityEffect(const EffectInfo * pEffectInfo);

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

#endif /*INOUTDEGREEASSORTATIVITYEFFECT_H_*/
