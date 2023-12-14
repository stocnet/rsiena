/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorChange.h
 *
 * Description: This file contains the definition of the BehaviorChange class.
 *****************************************************************************/


#ifndef BEHAVIORCHANGE_H_
#define BEHAVIORCHANGE_H_

#include <vector>

#include "MiniStep.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a ministep changing a behavior variable.
 */
class BehaviorChange: public MiniStep
{
public:
	BehaviorChange(BehaviorLongitudinalData * pData,
		int ego,
		int difference);
	virtual ~BehaviorChange();

	inline int difference() const;

	virtual bool behaviorMiniStep() const;
	virtual void makeChange(DependentVariable * pVariable);
	virtual bool missing(int period) const;
	virtual bool missingStart(int period) const;
	virtual bool missingEnd(int period) const;
	virtual MiniStep * createReverseMiniStep() const;
	virtual MiniStep * createCopyMiniStep() const;
	virtual bool firstOfConsecutiveCancelingPair() const;

private:
	// The longitudinal data object for the corresponding behavior variable
	BehaviorLongitudinalData * lpData;

	// The amount of change
	int ldifference{};

};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the amount of change in this ministep.
 */
int BehaviorChange::difference() const
{
	return this->ldifference;
}


}

#endif /* BEHAVIORCHANGE_H_ */
