/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkConstraint.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkConstraint.
 *****************************************************************************/

#include "NetworkConstraint.h"

namespace siena
{

/**
 * Constructs a new network constraint.
 */
NetworkConstraint::NetworkConstraint(std::string networkName1,
	std::string networkName2, NetworkConstraintType type)
{
	this->lnetworkName1 = networkName1;
	this->lnetworkName2 = networkName2;
	this->ltype = type;
}

}
