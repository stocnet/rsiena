/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoPathTable.h
 *
 * Description: This file defines the class TwoPathTable.
 *****************************************************************************/

#ifndef TWOPATHTABLE_H_
#define TWOPATHTABLE_H_

#include "EgocentricConfigurationTable.h"
#include "network/NetworkUtils.h"

namespace siena
{

/**
 * This class defines a table of generalized two-paths between a fixed actor,
 * the ego, and all other actors. In a generalized two-path, ties can be
 * traversed in either direction, or one can require to traverse only
 * reciprocated ties. For instance, one can think of in-stars
 * between two actors i and j as generalized two-paths, where we traverse
 * one outgoing tie of i, say (i,h), in the forward direction followed by the
 * incoming tie (j,h) of h in the backward direction.
 */
class TwoPathTable : public EgocentricConfigurationTable
{
public:
	TwoPathTable(NetworkCache * pOwner,
		Direction firstStepDirection,
		Direction secondStepDirection);

protected:
	virtual void calculate();

private:
	template<class Iterator> void performFirstStep(Iterator iter);
	template<class Iterator> void performSecondStep(Iterator iter);

	// The direction of the first step
	Direction lfirstStepDirection {};

	// The direction of the second step
	Direction lsecondStepDirection {};
};

}

#endif /*TWOPATHTABLE_H_*/
