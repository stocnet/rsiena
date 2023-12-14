/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkChange.h
 *
 * Description: This file contains the definition of the NetworkChange class.
 *****************************************************************************/

#ifndef NETWORKCHANGE_H_
#define NETWORKCHANGE_H_

#include <vector>

#include "MiniStep.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a ministep changing a network variable.
 */
class NetworkChange: public MiniStep
{
public:
	NetworkChange(NetworkLongitudinalData * pData,
		int ego,
		int alter,
		bool diagonal);
	virtual ~NetworkChange();

	virtual bool networkMiniStep() const;
	inline int alter() const;
	virtual void makeChange(DependentVariable * pVariable);
	virtual bool missing(int period) const;
	virtual bool missingStart(int period) const;
	virtual bool missingEnd(int period) const;
	virtual MiniStep * createReverseMiniStep() const;
	virtual MiniStep * createCopyMiniStep() const;

private:
	// The longitudinal data object for the corresponding network variable
	NetworkLongitudinalData * lpData;

	// The alter whose incoming tie is changed
	int lalter {};

};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the alter of this ministep.
 */
int NetworkChange::alter() const
{
	return this->lalter;
}

}

#endif /* NETWORKCHANGE_H_ */
