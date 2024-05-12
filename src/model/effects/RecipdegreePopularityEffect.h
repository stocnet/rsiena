/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecipdegreePopularityEffect.h
 *
 * Description: This file contains the definition of the
 * RecipdegreePopularityEffect class.
 *****************************************************************************/

#ifndef RECIPDEGREEPOPULARITYEFFECT_H_
#define RECIPDEGREEPOPULARITYEFFECT_H_

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
 * This class defines the reciprocal degree popularity effect defined by
 * s_i(x)= sum_j x_{ij} x^{(r)}_j, where
 * x^{(r)}_j is the reciprocated degree.
 */
class RecipdegreePopularityEffect : public NetworkEffect
{
public:
	RecipdegreePopularityEffect(const EffectInfo * pEffectInfo, bool root);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of indegrees must be used
	bool lroot {};

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*INDEGREEPOPULARITYEFFECT_H_*/
