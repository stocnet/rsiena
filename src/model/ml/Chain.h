/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.h
 *
 * Description: This file contains the definition of the Chain class.
 *****************************************************************************/

#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include <map>
#include "model/ml/Option.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Data;
class State;
class MLSimulation;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Defines a sequence of ministeps, which is the basic structure for
 * Maximum-Likelihood calculations.
 */
class Chain
{
public:
	Chain(Data * pData);
	virtual ~Chain();

	// Initialization

	void setupInitialState(bool copyValues);

	// Chain modifications

	void clear();
	void insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep);
	void remove(MiniStep * pMiniStep);
	void connect(int period, MLSimulation * pMLSimulation);
	void onReciprocalRateChange(const MiniStep * pMiniStep, double newValue);
	void changeInitialState(const MiniStep * pMiniStep);
	void recreateInitialState();
	void createInitialStateDifferences();
	void addInitialStateDifference(MiniStep * pMiniStep);
	void addEndStateDifference(MiniStep * pMiniStep);
	void clearEndStateDifferences();

	// Accessors

	void period(int period);
	int period() const;
	MiniStep * pFirst() const;
	MiniStep * pLast() const;
	const State * pInitialState() const;
	int ministepCount() const;
	int diagonalMinistepCount() const;
	int consecutiveCancelingPairCount() const;
	int missingNetworkMiniStepCount() const;
	int missingBehaviorMiniStepCount() const;
	double mu() const;
	double sigma2() const;
	double finalReciprocalRate() const;
	void finalReciprocalRate(double value);
	const std::vector <MiniStep *> & rEndStateDifferences() const;
	const std::vector <MiniStep *> & rInitialStateDifferences() const;

	// Intervals

	int intervalLength(const MiniStep * pFirstMiniStep,
		const MiniStep * pLastMiniStep) const;

	// Same option related

	MiniStep * firstMiniStepForOption(const Option & rOption) const;
	MiniStep * nextMiniStepForOption(const Option & rOption,
		const MiniStep * pFirstMiniStep) const;

	// Same dyad related

	MiniStep * pFirstMiniStepForLink(const MiniStep * pLinkMiniStep) const;
	MiniStep * pLastMiniStepForLink(const MiniStep * pLinkMiniStep) const;

	// Random draws

	MiniStep * randomMiniStep() const;
	MiniStep * randomDiagonalMiniStep() const;
	MiniStep * randomMiniStep(MiniStep * pFirstMiniStep,
		MiniStep * pLastMiniStep) const;
	MiniStep * randomConsecutiveCancelingPair() const;
	MiniStep * randomMissingNetworkMiniStep() const;
	MiniStep * randomMissingBehaviorMiniStep() const;

	// Copy
	Chain * copyChain() const;
//	void dumpChain() const;

private:
	void resetOrderingKeys();
	void updateSameOptionPointersOnInsert(MiniStep * pMiniStep);
	void updateCCPs(MiniStep * pMiniStep);

	// A dummy first ministep in the chain
	MiniStep * lpFirst;

	// A dummy last ministep in the chain
	MiniStep * lpLast;

	// The underlying observed data
	Data * lpData;

	// The period of changes represented by this chain
	int lperiod {};

	// The initial state of the variables (denoted y_init in the specification)
	State * lpInitialState;

	// The initial state of the variables stored as a vector of ministeps
	// relative to data (to save space)
	std::vector<MiniStep *> linitialStateDifferences;

	// The final state of the variables stored as a vector of ministeps
	// relative to data (to save space)
	std::vector<MiniStep *> lendStateDifferences;

	// Stores the ministeps in no particular order.
	// The first (dummy) ministep is not stored in this vector.

	std::vector<MiniStep *> lminiSteps;

	// Stores pointers to the diagonal ministeps in no particular order.
	std::vector<MiniStep *> ldiagonalMiniSteps;

	// Stores pointers to the first ministep of each CCP in no particular order.
	std::vector<MiniStep *> lccpMiniSteps;

	// Stores pointers to the networks ministeps of missing options
	std::vector<MiniStep *> lmissingNetworkMiniSteps;

	// Stores pointers to the behavior ministeps of missing options
	std::vector<MiniStep *> lmissingBehaviorMiniSteps;

	// Sum of reciprocal rates over all non-dummy ministeps.
	double lmu {};

	// Sum of squared reciprocal rates over all non-dummy ministeps.
	double lsigma2 {};

	// Final reciprocal rate when updating part or all of the probabilities.
	double lfinalReciprocalRate {};

	// Maps each option to its first ministep in the chain (if any)
	std::map<const Option, MiniStep *> lfirstMiniStepPerOption;
};

}

#endif /* CHAIN_H_ */
