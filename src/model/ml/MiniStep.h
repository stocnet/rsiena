/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.h
 *
 * Description: This file contains the definition of the MiniStep class.
 *****************************************************************************/


#ifndef MINISTEP_H_
#define MINISTEP_H_

#include <string>
#include <vector>
#include <map>

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class LongitudinalData;
class DependentVariable;
class Chain;
class Option;
class EffectInfo;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a single ministep as part of a Chain object in the
 * Maximum-Likelihood calculations.
 */
class MiniStep
{
	friend class Chain;

public:
	MiniStep(LongitudinalData * pData, int ego);
	virtual ~MiniStep();

	inline int ego() const;
	int variableId() const;
	std::string variableName() const;

	virtual bool networkMiniStep() const;
	virtual bool behaviorMiniStep() const;

	inline const Option * pOption() const;

	Chain * pChain() const;

	inline double logOptionSetProbability() const;
	void logOptionSetProbability(double probability);

	inline double logChoiceProbability() const;
	void logChoiceProbability(double probability);

	inline double reciprocalRate() const;
	void reciprocalRate(double value);

	inline MiniStep * pPrevious() const;
	inline MiniStep * pNext() const;
	void pPrevious(MiniStep * pMiniStep);
	void pNext(MiniStep * pMiniStep);

	inline MiniStep * pPreviousWithSameOption() const;
	inline MiniStep * pNextWithSameOption() const;
	void pPreviousWithSameOption(MiniStep * pMiniStep);
	void pNextWithSameOption(MiniStep * pMiniStep);

	inline int index() const;
	inline void index(int index);

	inline int diagonalIndex() const;
	inline void diagonalIndex(int index);

	inline int consecutiveCancelingPairIndex() const;
	inline void consecutiveCancelingPairIndex(int index);

	inline int missingIndex() const;
	inline void missingIndex(int index);

	inline double orderingKey() const;
	inline void orderingKey(double key);

	virtual void makeChange(DependentVariable * pVariable);
    bool diagonal() const;
	virtual bool missing(int period) const;
	virtual bool missingStart(int period) const;
	virtual bool missingEnd(int period) const;
	virtual MiniStep * createReverseMiniStep() const;
	virtual MiniStep * createCopyMiniStep() const;

	virtual bool firstOfConsecutiveCancelingPair() const;

	std::map<const EffectInfo *, std::vector<double> >* changeContributions() const;
	void changeContributions(std::map<const EffectInfo *, std::vector<double> > * contributions);


protected:
	void pOption(const Option * pOption);
    void diagonal(bool value);

private:
	void pChain(Chain * pChain);

	// The actor making the change
	int lego {};

	// The longitudinal data object for the corresponding dependent variable
	LongitudinalData * lpData;

	// The option of this ministep
	const Option * lpOption;

	// The owner chain (0, if the ministep is not part of a chain)
	Chain * lpChain;

	// Log probability of choosing the option set of this ministep,
	// given the state just before this ministep.

	double llogOptionSetProbability {};

	// Log probability of making this ministep, given that a ministep
	// of the same option set will be made.

	double llogChoiceProbability {};

	// Reciprocal of aggregate (summed) rate function immediately
	// before this ministep.

	double lreciprocalRate {};

	// Does this step result in a change of the dependent variables. Called
	// diagonal because originally diagonal entries in the adjacency matrices
	// were used to indicate such steps. Can occur with any actor and alter in
	// symmetric networks. In bipartite networks indicated by alter=number
	// of alters.
	bool ldiagonal {};

	// Points to the previous ministep in the same chain
	MiniStep * lpPrevious;

	// Points to the next ministep in the same chain
	MiniStep * lpNext;

	// Points to the previous ministep with the same option
	MiniStep * lpPreviousWithSameOption;

	// Points to the next ministep with the same option
	MiniStep * lpNextWithSameOption;

	// The index in the vector of ministeps of the owner chain.
	int lindex {};

	// The index in the vector of diagonal ministeps of the owner chain.
	int ldiagonalIndex {};

	// The index in the vector of CCPs of the owner chain.
	int lconsecutiveCancelingPairIndex {};

	// The index in the vector of network of behavior ministeps with
	// missing observed data at either end of the period

	int lmissingIndex {};

	// An arbitrary double value that corresponds with the chain order,
	// i.e. x.lorderingKey < y.lorderingKey whenever x and y is part
	// of the same chain and x precedes y in the chain.

	double lorderingKey {};

	// Stores for each effect its contributions to the tie flip probabilities or behavior change probabilities
	std::map<const EffectInfo *, std::vector<double> > * lpChangeContributions;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the ego of this ministep.
 */
int MiniStep::ego() const
{
	return this->lego;
}


/**
 * Returns the log probability of choosing the option set of this ministep,
 * given the state just before this ministep.
 */
double MiniStep::logOptionSetProbability() const
{
	return this->llogOptionSetProbability;
}


/**
 * Returns the option of this ministep.
 */
const Option * MiniStep::pOption() const
{
	return this->lpOption;
}


/**
 * Returns the log probability of making this ministep,
 * given that a ministep of the same option set will be made.
 */
double MiniStep::logChoiceProbability() const
{
	return this->llogChoiceProbability;
}


/**
 * Returns the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
double MiniStep::reciprocalRate() const
{
	return this->lreciprocalRate;
}


/**
 * Returns the previous ministep in the same chain.
 */
MiniStep * MiniStep::pPrevious() const
{
	return this->lpPrevious;
}


/**
 * Returns the next ministep in the same chain.
 */
MiniStep * MiniStep::pNext() const
{
	return this->lpNext;
}


/**
 * Returns the previous ministep having the same option as this ministep.
 */
MiniStep * MiniStep::pPreviousWithSameOption() const
{
	return this->lpPreviousWithSameOption;
}


/**
 * Returns the next ministep having the same option as this ministep.
 */
MiniStep * MiniStep::pNextWithSameOption() const
{
	return this->lpNextWithSameOption;
}


/**
 * Returns the index of this ministep in the vector of ministeps of the
 * owner chain.
 */
int MiniStep::index() const
{
	return this->lindex;
}


/**
 * Stores the index of this ministep in the vector of ministeps of the
 * owner chain.
 */
void MiniStep::index(int index)
{
	this->lindex = index;
}


/**
 * Returns the index of this ministep in the vector of diagonal
 * ministeps of the owner chain (-1, if the ministep is not part of a
 * chain or is not diagonal).
 */
int MiniStep::diagonalIndex() const
{
	return this->ldiagonalIndex;
}


/**
 * Stores the index of this ministep in the vector of diagonal
 * ministeps of the owner chain (-1, if the ministep is not part of a
 * chain or is not diagonal).
 */
void MiniStep::diagonalIndex(int index)
{
	this->ldiagonalIndex = index;
}


/**
 * Returns the index of this ministep in the vector of ministeps
 * of the owner chain that are the first ministeps of CCPs
 * (-1, if the ministep is not part of a chain or is not the first
 * ministep of a CCP).
 */
int MiniStep::consecutiveCancelingPairIndex() const
{
	return this->lconsecutiveCancelingPairIndex;
}


/**
 * Stores the index of this ministep in the vector of ministeps
 * of the owner chain that are the first ministeps of CCPs
 * (-1, if the ministep is not part of a chain or is not the first
 * ministep of a CCP).
 */
void MiniStep::consecutiveCancelingPairIndex(int index)
{
	this->lconsecutiveCancelingPairIndex = index;
}


/**
 * Returns the index of this ministep in the corresponding vector of
 * ministeps with missing observed data.
 * (-1, if the observed data is not missing).
 */
int MiniStep::missingIndex() const
{
	return this->lmissingIndex;
}


/**
 * Stores the index of this ministep in the corresponding vector of
 * ministeps with missing observed data.
 * (-1, if the observed data is not missing).
 */
void MiniStep::missingIndex(int index)
{
	this->lmissingIndex = index;
}

/**
 * Returns the ordering key of this ministep.
 */
double MiniStep::orderingKey() const
{
	return this->lorderingKey;
}


/**
 * Stores the ordering key of this ministep.
 */
void MiniStep::orderingKey(double key)
{
	this->lorderingKey = key;
}

}

#endif /* MINISTEP_H_ */
