/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DisjointFilter.h
 *
 * Description: This file contains the definition of the
 * DisjointFilter class.
 *****************************************************************************/

#ifndef DISJOINTFILTER_H_
#define DISJOINTFILTER_H_

#include "NetworkDependentFilter.h"

namespace siena
{

/**
 * Defines a filter of permissible changes that requires that every
 * tie (i,j) is present in no more than one of two networks.
 */
class DisjointFilter: public NetworkDependentFilter
{
public:
	DisjointFilter(const NetworkVariable * pOwnerVariable,
		const NetworkVariable * pOtherVariable);

	virtual void filterPermittedChanges(int ego, bool * permitted);
	virtual bool validMiniStep(const NetworkChange * pMiniStep);
	
private:
	bool lsymm {};
};

}

#endif /* DISJOINTFILTER_H_ */
