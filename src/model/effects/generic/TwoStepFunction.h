/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoStepFunction.h
 *
 * Description: This file contains the definition of the
 * TwoStepFunction class.
 *****************************************************************************/

#ifndef TWOSTEPFUNCTION_H_
#define TWOSTEPFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"
#include "network/NetworkUtils.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of two-paths between
 * the ego and alters in a network of the given name.
 */
class TwoStepFunction:
	public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	TwoStepFunction(std::string networkName, Direction direction1, Direction direction2);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter) const;
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;

	Direction ldirection1 {};
	Direction ldirection2 {};

};

}

#endif /* TWOSTEPFUNCTION_H_ */
