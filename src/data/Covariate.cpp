/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Covariate.cpp
 *
 * Description: This file contains the implementation of the
 * Covariate class.
 *****************************************************************************/

#include "Covariate.h"
#include "ActorSet.h"
#include <cmath>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction
// ----------------------------------------------------------------------------

/**
 * Creates a new covariate with the given name for the given set of actors.
 */
Covariate::Covariate(string name, const ActorSet * pActorSet) :
	NamedObject(name)
{
	this->lpActorSet = pActorSet;
	this->lmean = 0;
	this->lrange = 0;
	this->lsimilarityMean = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the domain set of actors for this covariate.
 */
const ActorSet * Covariate::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Stores the average value of this covariate.
 */
void Covariate::mean(double value)
{
	this->lmean = value;
}

/**
 * Stores the range of values of this covariate, which is calculated in R.
 */
void Covariate::range(double range)
{
	this->lrange = range;
}


int Covariate::covariateN() const
{
	return this->pActorSet()->n();
}

double Covariate::min() const
{
	double value = 0;
	return value;
}

double Covariate::max() const
{
	double value = 0;
	return value;
}

/**
 * Stores the similarity mean of this covariate, which is calculated in R.
 */
void Covariate::similarityMean(double similarityMean)
{
	this->lsimilarityMean = similarityMean;
}

/**
 * Stores the similarity alter mean of this covariate wrt this network,
 * which is calculated in R.
 */
void Covariate::similarityMeans(double similarityMean, string networkName)
{
	this->lsimilarityMeans[networkName] = similarityMean;
}

// ----------------------------------------------------------------------------
// Section: Similarity
// ----------------------------------------------------------------------------

/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double Covariate::similarity(double a, double b) const
{
	return 1.0 - fabs(a - b) / this->lrange - this->lsimilarityMean;
}

/**
 * Returns the centered alter similarity for the given values relative to the given
 * network, defined as 1 - |a - b| / range - similarityMeans[network].
 */
double Covariate::similarityNetwork(double a, double b, string networkName) const
{
	double similarityMean = 0;
	map<string, double>::const_iterator iter =
		this->lsimilarityMeans.find(networkName);
	if (iter != this->lsimilarityMeans.end())
	{
		similarityMean = iter->second;
	}
	return 1.0 - fabs(a - b) / this->lrange - similarityMean;
}

}
