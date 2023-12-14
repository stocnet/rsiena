/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkConstraint.h
 *
 * Description: This file contains the definition of the
 * NetworkConstraint class.
 *****************************************************************************/

#ifndef NETWORKCONSTRAINT_H_
#define NETWORKCONSTRAINT_H_

#include <string>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * This enumeration defines two possible directions of traversing a tie.
 */
enum NetworkConstraintType {HIGHER, DISJOINT, AT_LEAST_ONE};


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * A descriptor class for network constraints like higher(network1, network2).
 * These constraints are maintained during model simulations.
 */
class NetworkConstraint
{
public:
	NetworkConstraint(std::string networkName1,
		std::string networkName2,
		NetworkConstraintType type);

	inline std::string networkName1() const;
	inline std::string networkName2() const;
	inline NetworkConstraintType type() const;

private:
	// The name of the first network involved in this constraint
	std::string lnetworkName1 {};

	// The name of the second network involved in this constraint
	std::string lnetworkName2 {};

	// The type of the constraint
	NetworkConstraintType ltype {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the name of the first network variable involved in this constraint.
 */
std::string NetworkConstraint::networkName1() const
{
	return this->lnetworkName1;
}


/**
 * Returns the name of the second network variable involved in this constraint.
 */
std::string NetworkConstraint::networkName2() const
{
	return this->lnetworkName2;
}


/**
 * Returns the type of this constraint.
 */
NetworkConstraintType NetworkConstraint::type() const
{
	return this->ltype;
}

}

#endif /* NETWORKCONSTRAINT_H_ */
