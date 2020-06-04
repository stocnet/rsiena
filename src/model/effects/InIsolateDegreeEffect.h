/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InIsolateDegreeEffect.h
 *
 * Description: This file contains the definition of the
 * InIsolateDegreeEffect class.
 *****************************************************************************/

#ifndef INISOLATEDEGREEEFFECT_H_
#define INISOLATEDEGREEEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the in-isolate Outdegree effect defined by
 * s_i(x)= I{x_{+i}=0} x_{i+} 
 * The corresponding statistic is
 * the sum of outdegrees over all actors with indegree zero.
 */
class InIsolateDegreeEffect : public NetworkEffect
{
public:
	InIsolateDegreeEffect(const EffectInfo * pEffectInfo);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

};

}

#endif /*INISOLATEDEGREEEFFECT_H_*/
