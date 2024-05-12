/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConstantCovariate.h
 * 
 * Description: This file contains the definition of the ConstantCovariate
 * class.
 *****************************************************************************/

#ifndef CONSTANTCOVARIATE_H_
#define CONSTANTCOVARIATE_H_

#include "data/Covariate.h"

namespace siena
{

/**
 * This class stores the values of a covariate, which remains constant
 * throughout all observations.
 */
class ConstantCovariate : public Covariate
{
public:
	ConstantCovariate(std::string name, const ActorSet * pActorSet);
	virtual ~ConstantCovariate();
	
	double value(int i) const;
	void value(int i, double value);
	bool missing(int i) const;
	void missing(int i, bool missing);
	double min() const;
	double max() const;

private:
	// The values of this covariate
	double * lvalues {};
	double lmin {};
	double lmax {};

	// Missingness indicators
	bool * lmissing {};
};

}

#endif /*CONSTANTCOVARIATE_H_*/
