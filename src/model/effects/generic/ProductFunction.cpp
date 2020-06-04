#include "ProductFunction.h"

namespace siena
{

ProductFunction::ProductFunction(AlterFunction * pFirstFunction,
	AlterFunction * pSecondFunction)
{
	this->lpFirstFunction = pFirstFunction;
	this->lpSecondFunction = pSecondFunction;
}


/**
 * Deallocates this function.
 */
ProductFunction::~ProductFunction()
{
	delete this->lpFirstFunction;
	delete this->lpSecondFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ProductFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstFunction->initialize(pData, pState, period, pCache);
	this->lpSecondFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling ProductFunction::value(...).
 */
void ProductFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFirstFunction->preprocessEgo(ego);
	this->lpSecondFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ProductFunction::value(int alter)
{
	double value = this->lpFirstFunction->value(alter);

	// Hopefully the following if-statement improves the efficiency as
	// the calculation of the second function can be potentially expensive.

	if (value)
	{
		value *= this->lpSecondFunction->value(alter);
	}

	return value;
}

}
