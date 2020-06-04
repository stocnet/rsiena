#include "AlterPredicate.h"

namespace siena
{

/**
 * Creates a new predicate.
 */
AlterPredicate::AlterPredicate()
{
	this->lego = -1;
	this->lperiod = -1;
}


/**
 * Initializes this predicate.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AlterPredicate::initialize(const Data * pData,
	State * pState, int period, Cache * pCache)
{
	this->lperiod = period;
}


/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling AlterPredicate::value(...).
 */
void AlterPredicate::preprocessEgo(int ego)
{
	this->lego = ego;
}

}
