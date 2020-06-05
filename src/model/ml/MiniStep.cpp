/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.cpp
 *
 * Description: This file contains the implementation of the class MiniStep.
 *****************************************************************************/

#include "MiniStep.h"
#include "data/LongitudinalData.h"
#include "model/variables/DependentVariable.h"
#include "model/ml/Chain.h"

using namespace std;

namespace siena
{

/**
 * Constructs a new ministep.
 * @param[in] pData the longitudinal data object for the
 * corresponding dependent variable
 * @param[in] ego the actor making the change
 */
MiniStep::MiniStep(LongitudinalData * pData, int ego)
{
	this->lego = ego;
	this->lpData = pData;
	this->lpOption = 0;
	this->lpChain = 0;
	this->llogOptionSetProbability = 0;
	this->llogChoiceProbability = 0;
	this->lreciprocalRate = 0;
	this->lpPrevious = 0;
	this->lpNext = 0;
	this->lpPreviousWithSameOption = 0;
	this->lpNextWithSameOption = 0;
	this->lindex = -1;
	this->ldiagonalIndex = -1;
	this->lconsecutiveCancelingPairIndex = -1;
	this->lmissingIndex = -1;
	this->lorderingKey = 0;
	this->ldiagonal = false;
	this->lpChangeContributions = 0;
}


/**
 * Deallocates this ministep.
 */
MiniStep::~MiniStep()
{
	if (this->lpOption)
	{
		delete this->lpOption;
	}

	this->lpOption = 0;
	if (this->lpChangeContributions)
	{
		delete this->lpChangeContributions;
	}
}


/**
 * Returns if this ministep is changing a network variable.
 */
bool MiniStep::networkMiniStep() const
{
	return false;
}


/**
 * Returns if this ministep is changing a behavior variable.
 */
bool MiniStep::behaviorMiniStep() const
{
	return false;
}


/**
 * Returns the ID of the dependent variable that this ministep is changing.
 */
int MiniStep::variableId() const
{
	return this->lpData->id();
}

/**
 * Returns the name of the dependent variable that this ministep is changing.
 */
string MiniStep::variableName() const
{
	return this->lpData->name();
}




/**
 * Stores the owner chain of this ministep.
 */
void MiniStep::pChain(Chain * pChain)
{
	this->lpChain = pChain;
}


/**
 * Returns the owner chain of this ministep.
 */
Chain * MiniStep::pChain() const
{
	return this->lpChain;
}


/**
 * Stores the option of this ministep.
 */
void MiniStep::pOption(const Option * pOption)
{
	this->lpOption = pOption;
}

/**
 * Stores whether this ministep is diagonal.
 */
void MiniStep::diagonal(bool value)
{
	this->ldiagonal = value;
}

/**
 * Stores the log probability of choosing the option set of this ministep,
 * given the state just before this ministep.
 */
void MiniStep::logOptionSetProbability(double probability)
{
	this->llogOptionSetProbability = probability;
}


/**
 * Stores the log probability of making this ministep,
 * given that a ministep of the same option set will be made.
 */
void MiniStep::logChoiceProbability(double probability)
{
	this->llogChoiceProbability = probability;
}


/**
 * Stores the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
void MiniStep::reciprocalRate(double value)
{
	if (this->lpChain)
	{
		this->lpChain->onReciprocalRateChange(this, value);
	}

	this->lreciprocalRate = value;
}


/**
 * Stores the pointer to the previous ministep in the same chain.
 */
void MiniStep::pPrevious(MiniStep * pMiniStep)
{
	this->lpPrevious = pMiniStep;
}


/**
 * Stores the pointer to the next ministep in the same chain.
 */
void MiniStep::pNext(MiniStep * pMiniStep)
{
	this->lpNext = pMiniStep;
}


/**
 * Stores the pointer to the previous ministep having the same option as this
 * ministep.
 */
void MiniStep::pPreviousWithSameOption(MiniStep * pMiniStep)
{
	this->lpPreviousWithSameOption = pMiniStep;
}


/**
 * Stores the pointer to the next ministep having the same option as this
 * ministep.
 */
void MiniStep::pNextWithSameOption(MiniStep * pMiniStep)
{
	this->lpNextWithSameOption = pMiniStep;
}


/**
 * Changes the given dependent variable according to this ministep.
 */
void MiniStep::makeChange(DependentVariable * pVariable)
{
	// Nothing in the base class.
}


/**
 * Returns if this ministep is diagonal, namely, it does not change
 * the dependent variables. Dummy ministeps are not considered diagonal.
 */
bool MiniStep::diagonal() const
{
	return this->ldiagonal;
}


/**
 * Returns if the observed data for this ministep is missing at
 * either end of the given period.
 */
bool MiniStep::missing(int period) const
{
	return false;
}
/**
 * Returns if the observed data for this ministep is missing at
 * the start of the given period.
 */
bool MiniStep::missingStart(int period) const
{
	return false;
}

/**
 * Returns if the observed data for this ministep is missing at
 * the end of the given period.
 */
bool MiniStep::missingEnd(int period) const
{
	return false;
}

/**
 * Returns a new ministep that reverses the effect of this ministep.
 */
MiniStep * MiniStep::createReverseMiniStep() const
{
	return 0;
}

/**
 * Returns a new ministep that is a copy of this ministep.
 */
MiniStep * MiniStep::createCopyMiniStep() const
{
	return 0;
}

/**
 * Returns if this mini step is the first mini step of a CCP.
 */
bool MiniStep::firstOfConsecutiveCancelingPair() const
{
	//Rprintf(" first ccp %x\n", this->pChain());
	bool missing = false;
	if (this->pChain())
	{
		missing = this->missing(this->pChain()->period());
	}
	return
		!this->diagonal() &&
		this->lpNextWithSameOption &&
		this->lpNextWithSameOption != this->lpNext &&
		!missing;
}

/**
 * Stores the contributions of each effect to possible
 * tie flips or behavior changes before this ministep.
 */
void MiniStep::changeContributions(map<const EffectInfo *, vector<double> > * contributions)
{
	this->lpChangeContributions = contributions;
}

/**
 * Returns the contributions of each effect to possible
 * tie flips or behavior changes before this ministep.
 */
map<const EffectInfo *, vector<double> > * MiniStep::changeContributions() const
{
	return this->lpChangeContributions;
}

}
