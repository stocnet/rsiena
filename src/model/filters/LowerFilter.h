/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LowerFilter.h
 *
 * Description: This file contains the definition of the LowerFilter class.
 *****************************************************************************/

#ifndef LOWERFILTER_H_
#define LOWERFILTER_H_

#include "NetworkDependentFilter.h"

namespace siena
{

/**
 * Defines a filter of permissible changes that requires that
 * x_{ij} <= y_{ij} for every all pairs (i,j) and two networks X and Y.
 */
class LowerFilter: public NetworkDependentFilter
{
public:
	LowerFilter(const NetworkVariable * pOwnerVariable,
		const NetworkVariable * pOtherVariable);

	virtual void filterPermittedChanges(int ego, bool * permitted);
	virtual bool validMiniStep(const NetworkChange * pMiniStep);
};

}

#endif /* LOWERFILTER_H_ */
