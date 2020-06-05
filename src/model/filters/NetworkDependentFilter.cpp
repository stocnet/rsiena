/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentFilter.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkDependentFilter.
 *****************************************************************************/

#include "NetworkDependentFilter.h"

namespace siena
{

/**
 * Constructs a new filter.
 */
NetworkDependentFilter::NetworkDependentFilter(
	const NetworkVariable * pOwnerVariable,
	const NetworkVariable * pOtherVariable) :
		PermittedChangeFilter(pOwnerVariable)
{
	this->lpOtherVariable = pOtherVariable;
}

}
