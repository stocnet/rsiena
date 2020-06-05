/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentFilter.h
 *
 * Description: This file contains the definition of the
 * NetworkDependentFilter class.
 *****************************************************************************/

#ifndef NETWORKDEPENDENTFILTER_H_
#define NETWORKDEPENDENTFILTER_H_

#include "PermittedChangeFilter.h"

namespace siena
{

/**
 * Defines a filter of permissible changes that depends on another network.
 */
class NetworkDependentFilter: public PermittedChangeFilter
{
public:
	NetworkDependentFilter(const NetworkVariable * pOwnerVariable,
		const NetworkVariable * pOtherVariable);

protected:
	inline const NetworkVariable * pOtherVariable() const;

private:
	// The other network that this filter depends on
	const NetworkVariable * lpOtherVariable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the (other) network variable that this filter depends on.
 */
const NetworkVariable * NetworkDependentFilter::pOtherVariable() const
{
	return this->lpOtherVariable;
}

}

#endif /* NETWORKDEPENDENTFILTER_H_ */
