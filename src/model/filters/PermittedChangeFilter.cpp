/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PermittedChangeFilter.cpp
 *
 * Description: This file contains the implementation of the class
 * PermittedChangeFilter.
 *****************************************************************************/

#include "PermittedChangeFilter.h"

namespace siena
{

/**
 * Creates a new filter of permissible changes for the given variable.
 */
PermittedChangeFilter::PermittedChangeFilter(const NetworkVariable * pVariable)
{
	this->lpVariable = pVariable;
}

}
