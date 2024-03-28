/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ChangingCovariate.h
 * 
 * Description: This file defines the ChangingCovariate class.
 *****************************************************************************/

#ifndef CHANGINGCOVARIATE_H_
#define CHANGINGCOVARIATE_H_

#include "data/Covariate.h"

namespace siena
{

/**
 * This class stores the values of a covariate, which is changing from
 * one observation to another.
 */
class ChangingCovariate : public Covariate
{
public:
	ChangingCovariate(std::string name,
		const ActorSet * pActorSet,
		int observationCount);
	virtual ~ChangingCovariate();
	
	double value(int i, int observation) const;
	void value(int i, int observation, double value);
	bool missing(int i, int observation) const;
	void missing(int i, int observation, bool missing);
	double min() const;
	double max() const;

private:
	// The values of the covariate per each actor and observation
	double ** lvalues {};
	double lmin {};
	double lmax {};
	
	// Missingness indicators
	bool ** lmissing {};
};

}

#endif /*CHANGINGCOVARIATE_H_*/
