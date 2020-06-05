/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AtLeastOneFilter.h
 *
 * Description: This file contains the definition of the
 * AtLeastOneFilter class.
 *****************************************************************************/

#ifndef ATLEASTONEFILTER_H_
#define ATLEASTONEFILTER_H_

#include "NetworkDependentFilter.h"

namespace siena
{

/**
 * Defines a filter of permissible changes that requires that every
 * tie (i,j) is present in at least one of two networks.
 */
class AtLeastOneFilter: public NetworkDependentFilter
{
public:
	AtLeastOneFilter(const NetworkVariable * pOwnerVariable,
		const NetworkVariable * pOtherVariable);

	virtual void filterPermittedChanges(int ego, bool * permitted);
	virtual bool validMiniStep(const NetworkChange * pMiniStep);
};

}

#endif /* ATLEASTONEFILTER_H_ */
