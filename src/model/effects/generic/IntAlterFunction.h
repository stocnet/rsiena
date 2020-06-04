/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntAlterFunction.h
 *
 * Description: This file contains the definition of the
 * IntAlterFunction class.
 *****************************************************************************/


#ifndef INTALTERFUNCTION_H_
#define INTALTERFUNCTION_H_

namespace siena
{

class Data;
class State;
class Cache;


/**
 * An interface for alter functions with integer values.
 */
class IntAlterFunction
{
public:
	/**
	 * Returns the value of this function for the given alter.
	 */
	virtual int intValue(int alter) = 0;
};

}

#endif /* INTALTERFUNCTION_H_ */
