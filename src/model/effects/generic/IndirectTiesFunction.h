/******************************************************************************
 * SIENA: Simulation Outvestigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File:IndirectTiesFunction.h
 *
 * Description: This file contains the definition of the
 *IndirectTiesFunction class.
 *****************************************************************************/


#ifndef INDIRECTTIESFUNCTION_H_
#define INDIRECTTIESFUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

class ConfigurationTable;

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;


// ----------------------------------------------------------------------------
// Section: IndirectTiesFunction class
// ----------------------------------------------------------------------------


/**
 * Defines a function that returns the number of actors at distance 2
 * from ego in a network of the given name.
 */
class IndirectTiesFunction: public NetworkAlterFunction
{
public:
	IndirectTiesFunction(std::string networkName, double parameter, 
				bool firstIn, bool secondIn);
	virtual ~IndirectTiesFunction();


	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);

	virtual double value(int alter) const; 

private:
	bool lroot {}; // should the square root be taken?
	bool lfirstin {}; // first tie in- or outgoing?
	bool lsecondin {}; // second tie in- or outgoing?
	// Lookup table for fast square root calculations:
	SqrtTable * lsqrtTable;
	// A helper array of marks for statistic calculation
	bool * lmark {};
	int lNbrDist2Nodes {};
};

}

#endif /* INDIRECTTIESFUNCTION_H_ */
