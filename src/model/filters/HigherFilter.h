/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HigherFilter.h
 *
 * Description: This file contains the definition of the HigherFilter class.
 *****************************************************************************/

#ifndef HIGHERFILTER_H_
#define HIGHERFILTER_H_

#include "NetworkDependentFilter.h"
#include <string>

namespace siena
{

/**
 * Defines a filter of permissible changes that requires that
 * x_{ij} >= y_{ij} for every all pairs (i,j) and two networks X and Y.
 */
class HigherFilter: public NetworkDependentFilter
{
public:
	HigherFilter(const NetworkVariable * pOwnerVariable,
		const NetworkVariable * pOtherVariable);

	virtual void filterPermittedChanges(int ego, bool * permitted);
	virtual bool validMiniStep(const NetworkChange * pMiniStep);
	
private:
	bool lsymm {};
};

}

#endif /* HIGHERFILTER_H_ */
