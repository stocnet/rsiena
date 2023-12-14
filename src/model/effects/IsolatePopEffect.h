/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IsolatePopEffect.h
 *
 * Description: This file contains the definition of the
 * IsolatePopEffect class.
 *****************************************************************************/

#ifndef ISOLATEPOPEFFECT_H_
#define ISOLATEPOPEFFECT_H_

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
 * This class defines the isolate popularity effect defined by
 * s_i(x)=  \sum_j x_{ij} I\{x_{+j} = 1, x_{j+} = 0 }
 * The corresponding statistic is
 * the number of ties to alters with indegree 1 and outdegree 0.
 */
class IsolatePopEffect : public NetworkEffect
{
public:
	IsolatePopEffect(const EffectInfo * pEffectInfo, bool outgoing);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// outgoing indicates whether the outdegree=0 of alter also
	// is taken into account.
	// outgoing=false was used for the shortlived inIsolatePop effect
	// may be dropped now
	bool loutgoing {};
};

}

#endif /*ISOLATEPOPEFFECT_H_*/
