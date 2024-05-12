/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AlterPredicate.h
 *
 * Description: This file contains the definition of the AlterPredicate class.
 *****************************************************************************/

#ifndef ALTERPREDICATE_H_
#define ALTERPREDICATE_H_

namespace siena {

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class State;
class Cache;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a predicate on all actors (alters) with respect to a given ego.
 */
class AlterPredicate {
public:
	virtual ~AlterPredicate() {}
	virtual void initialize(const Data * pData, State * pState, int period,
			Cache * pCache);
	virtual void preprocessEgo(int ego);

	inline int ego() const;
	inline int period() const;

	/**
	 * Returns if this predicate holds for the given alter. It is assumed
	 * that the predicate has been initialized before and pre-processed with
	 * respect to a certain ego.
	 */
	virtual bool value(int alter) = 0;

protected:
	AlterPredicate();

private:
	int lego{};
	int lperiod{};
};

// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the current ego.
 */
int AlterPredicate::ego() const {
	return this->lego;
}

/**
 * Returns the period of interest.
 */
int AlterPredicate::period() const {
	return this->lperiod;
}

}

#endif /* ALTERPREDICATE_H_ */
