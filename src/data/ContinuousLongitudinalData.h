/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ContinuousLongitudinalData.h
 *
 * Description: This module defines the class ContinuousLongitudinalData.
 *****************************************************************************/

#ifndef CONTINUOUSLONGITUDINALDATA_H_
#define CONTINUOUSLONGITUDINALDATA_H_

#include <map>

#include "data/LongitudinalData.h"

namespace siena
{

/**
 * This class stores the observed values of a continuous behavioral variable 
 * for one or more observation moments.
 */
class ContinuousLongitudinalData : public LongitudinalData
{
public:
	ContinuousLongitudinalData(int id,
		std::string name,
		const ActorSet * pActorSet,
		int observationCount);
	virtual ~ContinuousLongitudinalData();

	double value(int observation, int actor) const;
	void value(int observation, int actor, double value);
	bool missing(int observation, int actor) const;
	void missing(int observation, int actor, bool missing);
	bool structural(int observation, int actor) const;
	void structural(int observation, int actor, bool flag);
	const double * values(int observation) const;
	const double * valuesLessMissings(int observation) const;
	const double * valuesLessMissingStarts(int observation) const;

	double min() const;
	double max() const;
	double overallMean() const;
	double range() const;
	double similarity(double a, double b) const;
	double similarityNetwork(double a, double b, std::string networkName) const;
	double similarityMean() const;
	void similarityMean(double similarityMean);
	void similarityMeans(double similarityMean, std::string networkName);
	virtual double observedDistribution(int value, int observation) const;
	void calculateProperties();

private:
	// An array of values per each observation
	double ** lvalues {};

	// Missingness indicators
	bool ** lmissing {};

	// Structural value indicators
	bool ** lstructural {};

	// An array of values per each observation
	double ** lvaluesLessMissings {};

	// An array of values per each observation
	double ** lvaluesLessMissingStarts {};

	// The smallest non-missing value
	double lmin {};

	// The largest non-missing value
	double lmax {};

	// The overall mean value over all observations
	double loverallMean {};

	// The range of observed values
	double lrange {};

	// The similarity mean
	double lsimilarityMean {};

	// The alter similarity means for each network (to be passed from outside)
	std::map<std::string, double> lsimilarityMeans;
};

}

#endif /*CONTINUOUSLONGITUDINALDATA_H_*/
