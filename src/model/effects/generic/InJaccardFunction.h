/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InJaccardFunction.h
 *
 * Description: This file contains the definition of the
 * InJaccardFunction class.
 *****************************************************************************/


#ifndef INJACCARDFUNCTION_H_
#define INJACCARDFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of in-stars between
 * the ego and alters in a network of the given name.
 */
class InJaccardFunction: public NetworkAlterFunction
{
public:
	InJaccardFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

	virtual double value(int alter);

private:
	ConfigurationTable * lpTable;
};

}

#endif /* INJACCARDFUNCTION_H_ */
