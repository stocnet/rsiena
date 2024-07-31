/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AlterFunction.h
 *
 * Description: This file contains the definition of the
 * AlterFunction class.
 *****************************************************************************/


#ifndef ALTERFUNCTION_H_
#define ALTERFUNCTION_H_

namespace siena
{

class Data;
class State;
class Cache;

class AlterFunction
{
public:
	AlterFunction();
	virtual ~AlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void initialize(const Data * pData,
		State * pState, State * pSimulatedState, int period, Cache * pCache);
	virtual void preprocessEgo(int ego);

	inline int ego() const;

	/**
	 * Will return the value of this function for the given alter. It is assumed
	 * that the function has been initialized before and pre-processed with
	 * respect to a certain ego.
	 */
// the "virtual ... = 0" means that overriding this method is required
// a method constructed in this way is called a pure method.
	virtual double value(int alter) const = 0; 

private:
	int lego{};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the current ego.
 */
int AlterFunction::ego() const
{
	return this->lego;
}


}

#endif /* ALTERFUNCTION_H_ */
