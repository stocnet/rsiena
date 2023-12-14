/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorLongitudinalData.h
 *
 * Description: This module defines the class BehaviorLongitudinalData.
 *****************************************************************************/

#ifndef BEHAVIORLONGITUDINALDATA_H_
#define BEHAVIORLONGITUDINALDATA_H_

#include <map>

#include "data/LongitudinalData.h"

namespace siena
{

/**
 * This class stores the observed values of a behavior variable for one or more
 * observation moments.
 */
class BehaviorLongitudinalData : public LongitudinalData
{
public:
	BehaviorLongitudinalData(int id,
		std::string name,
		const ActorSet * pActorSet,
		int observationCount);
	virtual ~BehaviorLongitudinalData();

	int value(int observation, int actor) const;
	void value(int observation, int actor, int value);
	bool missing(int observation, int actor) const;
	void missing(int observation, int actor, bool missing);
	bool structural(int observation, int actor) const;
	void structural(int observation, int actor, bool flag);
	const int * values(int observation) const;
	const int * valuesLessMissings(int observation) const;
	const int * valuesLessMissingStarts(int observation) const;
	void behModelType(int type);
	int behModelType() const;

	int min() const;
	int max() const;
	double overallMean() const;
	int range() const;
	double similarity(double a, double b) const;
	double similarityNetwork(double a, double b, std::string networkName) const;
	double similarityMean() const;
	void similarityMean(double similarityMean);
	void similarityMeans(double similarityMean, std::string networkName);
	virtual double observedDistribution(int value, int observation) const;
	void calculateProperties();

private:
	// An array of values per each observation
	int ** lvalues {};

	// Missingness indicators
	bool ** lmissing {};

	// Structural value indicators
	bool ** lstructural {};

	// An array of values per each observation
	int ** lvaluesLessMissings {};

	// An array of values per each observation
	int ** lvaluesLessMissingStarts {};

	// The behavioral model type.
	int lbehModelType {};

	// The smallest non-missing value
	int lmin {};

	// The largest non-missing value
	int lmax {};

	// The overall mean value over all observations (zamean in Pascal)
	double loverallMean {};

	// The range of observed values
	int lrange {};

	// The similarity mean
	double lsimilarityMean {};

	// The alter similarity means for each network (to be passed from outside)
	std::map<std::string, double> lsimilarityMeans;

	// The distribution of observed values at each observation.
	// lobservedDistributions[observation][value] stores the frequency of
	// the given value at the given observation.

	std::map<int, double> * lobservedDistributions;
};

}

#endif /*BEHAVIORLONGITUDINALDATA_H_*/
