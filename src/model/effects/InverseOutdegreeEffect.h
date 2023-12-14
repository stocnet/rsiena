/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InverseOutdegreeEffect.h
 *
 * Description: This file contains the declaration of the class
 * InverseOutdegreeEffect.
 *****************************************************************************/

#ifndef INVERSEOUTDEGREEEFFECT_H_
#define INVERSEOUTDEGREEEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the inverse outdegree effect defined as
 * s_i = 1/(outdegree(i) + c), where c is a parameter.
 * See the manual for effect definitions.
 */
class InverseOutdegreeEffect : public NetworkEffect
{
public:
	InverseOutdegreeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:
	double lc {};
};

}

#endif /*INVERSEOUTDEGREEEFFECT_H_*/
