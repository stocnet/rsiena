/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorChange.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorChange.
 *****************************************************************************/

#include "BehaviorChange.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/variables/BehaviorVariable.h"
#include "model/ml/Option.h"

namespace siena
{

/**
 * Constructs a new behavior ministep.
 * @param[in] pData the longitudinal data object for the
 * @param[in] ego the actor making the change
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
BehaviorChange::BehaviorChange(BehaviorLongitudinalData * pData,
	int ego,
	int difference) : MiniStep(pData, ego)
{
	this->lpData = pData;
	this->ldifference = difference;
	this->pOption(new Option(pData->id(), ego));
	this->diagonal(difference == 0);
}


/**
 * Deallocates this ministep.
 */
BehaviorChange::~BehaviorChange()
{
}


/**
 * Returns if this ministep is changing a behavior variable.
 */
bool BehaviorChange::behaviorMiniStep() const
{
	return true;
}


/**
 * Changes the given behavior variable according to this ministep.
 */
void BehaviorChange::makeChange(DependentVariable * pVariable)
{
	MiniStep::makeChange(pVariable);

	if (this->difference() != 0)
	{
		BehaviorVariable * pBehaviorVariable =
			dynamic_cast<BehaviorVariable *>(pVariable);
		int oldValue = pBehaviorVariable->value(this->ego());
		pBehaviorVariable->value(this->ego(), oldValue + this->difference());
	}
}


// /**
//  * Returns if this ministep is diagonal, namely, it does not change
//  * the dependent variables.
//  */
// bool BehaviorChange::diagonal() const
// {
// 	return this->difference() == 0;
// }


/**
 * Returns if the observed data for this ministep is missing at
 * either end of the given period.
 */
bool BehaviorChange::missing(int period) const
{
//	return this->lpData->missing(period, this->ego()) ||
//		this->lpData->missing(period + 1, this->ego());
//	return this->lpData->missing(period, this->ego()) ||
//		this->lpData->missing(period + 1, this->ego());
	return this->missingStart(period) || this->missingEnd(period);
}

/**
 * Returns if the observed data for this ministep is missing at
 * the start of the given period.
 */
bool BehaviorChange::missingStart(int period) const
{
//	return this->lpData->missing(period, this->ego()) ||
//		this->lpData->missing(period + 1, this->ego());
	return this->lpData->missing(period, this->ego());
}
/**
 * Returns if the observed data for this ministep is missing at
 * the end of the given period.
 */
bool BehaviorChange::missingEnd(int period) const
{
//	return this->lpData->missing(period, this->ego()) ||
//		this->lpData->missing(period + 1, this->ego());
	return this->lpData->missing(period + 1, this->ego());
}

/**
 * Returns a new ministep that reverses the effect of this ministep.
 */
MiniStep * BehaviorChange::createReverseMiniStep() const
{
	return new BehaviorChange(this->lpData,
		this->ego(),
		-this->difference());
}

/**
 * Returns a new ministep that is a copy of this ministep.
 */
MiniStep * BehaviorChange::createCopyMiniStep() const
{
	BehaviorChange * pBehaviorChange=  new BehaviorChange(this->lpData,
		this->ego(), this->difference());

	return pBehaviorChange;
}

/**
 * Returns if this mini step is the first mini step of a CCP.
 */
bool BehaviorChange::firstOfConsecutiveCancelingPair() const
{
	bool rc = MiniStep::firstOfConsecutiveCancelingPair();

	if (rc)
	{
		BehaviorChange * pNextForOption =
			dynamic_cast<BehaviorChange *>(this->pNextWithSameOption());

		rc = this->ldifference == -pNextForOption->ldifference;
	}

	return rc;
}


}
