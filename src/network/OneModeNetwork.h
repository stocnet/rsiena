/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetwork.h
 *
 * Description: This file contains the definition of the
 * OneModeNetwork class.
 *****************************************************************************/

#ifndef ONEMODENETWORK_H_
#define ONEMODENETWORK_H_

#include "Network.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class CommonNeighborIterator;
class UnionNeighborIterator;


// ----------------------------------------------------------------------------
// Section: OneModeNetwork class
// ----------------------------------------------------------------------------

/**
 * This class defines a regular network where a single set of actors acts as
 * both senders and receivers of ties.
 */
class OneModeNetwork : public Network
{
public:
	OneModeNetwork(int n, bool loopsPermitted = false);
	OneModeNetwork(const OneModeNetwork & rNetwork);
	OneModeNetwork & operator=(const OneModeNetwork & rNetwork);
	virtual Network * clone() const;
	virtual ~OneModeNetwork();

	bool loopsPermitted() const;
	int reciprocalDegree(int i) const;
	bool symmetric() const;
	virtual void clear();
	virtual bool isOneMode() const;

	CommonNeighborIterator reciprocatedTies(int i) const;
	CommonNeighborIterator reciprocatedTies(int i,
		int lowerBound) const;
	UnionNeighborIterator eitherTies(int i) const;

	int twoPathCount(int i, int j) const;
	int truncatedTwoPathCount(int i, int j, int threshold = 2) const;
	bool noTwoPaths(int i, int j, int intermediateActorUpperBound) const;
	bool existsTwoPath(int i, int j) const;
	bool atMostKTwoPaths(int i, int j, int k, int & twoPathCount) const;
	void neighborCensus(int i, int j, int & n3, int & n4) const;

protected:
	virtual int changeTieValue(int i, int j, int v, ChangeType type);
	virtual void onTieWithdrawal(int i, int j);
	virtual void onTieIntroduction(int i, int j);
	virtual int maxTieCount() const;

private:
	// Indicates if loops are permitted in this network
	bool lloopsPermitted {};

	// The reciprocal degree of each actor
	int * lpReciprocalDegree {};
};

}

#endif /*ONEMODENETWORK_H_*/
