/******************************************************************************
 * SIENA: Simulation Outvestigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutJaccardFunction.h
 *
 * Description: This file contains the definition of the
 * OutJaccardFunction class.
 *****************************************************************************/


#ifndef OUTJACCARDFUNCTION_H_
#define OUTJACCARDFUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of in-stars between
 * the ego and alters in a network of the given name.
 */
class OutJaccardFunction: public NetworkAlterFunction
{
public:
	OutJaccardFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);

private:
	ConfigurationTable * lpTable;
};

}

#endif /* OUTJACCARDFUNCTION_H_ */
