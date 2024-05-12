/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecipdegreeActivityEffect.h
 *
 * Description: This file contains the definition of the
 * RecipdegreeActivityEffect class.
 *****************************************************************************/

#ifndef RECIPDEGREEACTIVITYEFFECT_H_
#define RECIPDEGREEACTIVITYEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

/**
 * This class defines the reciprocal degree activity effect defined by
 * s_i(x)= sum_j x_{ij} x^{(r)}_i, where 
 * x^{(r)}_i is the reciprocated degree.
 */
class RecipdegreeActivityEffect : public NetworkEffect
{
public:
	RecipdegreeActivityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of reciprocated degrees must be used
	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*RECIPDEGREEACTIVITYEFFECT_H_*/
