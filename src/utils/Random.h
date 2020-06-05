/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Random.h
 *
 * Description: This module defines utilities for random number
 * generation.
 *****************************************************************************/

#ifndef RANDOM_H_
#define RANDOM_H_

#include <vector>

namespace siena
{

// Drawing random number from various distributions

double nextDouble();
double nextExponential(double lambda);
double nextExponentialQAD(double lambda);
double nextGamma(double shape, double scale);
double nextNormal(double mean, double standardDeviation);
int nextInt(int n);
int nextIntWithProbabilities(int n, const double * p);

// Density function for normal variables
double normalDensity(double value, double mean, double standardDeviation,
	int log);

//void storeState();
//void restoreState();


/**
 * Draws a random integer <i>i</i> from the interval [0,<i>n</i>)
 * according to the cumulative probability vector proportional to <i>p</i>,
 * namely, P(<i>i</i> <= <i>k</i>) = <i>p</i>[<i>k</i>]/<i>p</i>[<i>n</i>-1].
 */
template<class T>
int nextIntWithCumulativeProbabilities(int n, const T * p)
{
	// Take a uniform random number 'value' in (0, p[n - 1]) and perform
	// a binary search for i such that p[i] >= value and p[j] < value for
	// all j < i.

	double value = nextDouble() * p[n - 1];
	int i = 0;
	int high = n - 1;

	while (p[i] < value)
	{
		int middle = (i + high) / 2;

		if (p[middle] < value)
		{
			i = middle + 1;
		}
		else
		{
			high = middle;
		}
	}

	return i;
}


/**
 * Permutes the elements in the given vector.
 */
template<class T>
void permuteVector(std::vector<T> & rVector)
{
	for (unsigned i = 1; i < rVector.size(); i++)
	{
		T element = rVector[i];
		int newIndex = nextInt(i + 1);
		rVector[i] = rVector[newIndex];
		rVector[newIndex] = element;
	}
}

}

#endif /*RANDOM_H_*/
