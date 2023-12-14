#ifndef MLSIMULATION_H_
#define MLSIMULATION_H_

#include "model/EpochSimulation.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

enum Aspect {NETWORK, BEHAVIOR};

/**
 * This enumeration defines the possible types of MH step
 */

enum MHStepType {INSDIAG, CANCDIAG, PERMUTE, INSPERM, DELPERM, INSMISS,
				 DELMISS, INSMISDAT, DELMISDAT, MOVE, NBRTYPES};
// This is an implicit way of defining NBRTYPES as the number of step types, 10.
// INSMISDAT and DELMISDAT are versions of INSPERM and DELPERM
// that do not have their own probability for selecting a stepType.

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Option;
class DependentVariable;
class NetworkVariable;

/**
 * This class provides the functionality necessary for simulating an ML model
 * between two observations.
 */

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class MLSimulation: public EpochSimulation
{
public:
	MLSimulation(Data * pData, Model * pModel);
	virtual ~MLSimulation();

    void initialize(int period);
	void initializeInitialState(int period);
	void connect(int period);
	void preburnin();
	void runEpoch(int period);
	void MLStep();
	void setUpProbabilityArray();
    void updateProbabilities(Chain * pChain,
    	MiniStep * pFirstMiniStep,
    	MiniStep * pLastMiniStep);
    void executeMiniSteps(MiniStep * pFirstMiniStep, MiniStep * pLastMiniStep);

//	int acceptances(int stepType) const;
//	int rejections(int stepType) const;
	int aborts(int stepType) const;
	void incrementAborts(int stepType);

    // Metropolis-Hastings steps

	bool insertDiagonalMiniStep();
	bool cancelDiagonalMiniStep();
	bool permute(int c0);
	bool move();
	bool insertPermute(int c0);
	bool deletePermute(int c0);
	bool insertMissing();
	bool deleteMissing();
	double proposalProbability() const;
	bool missingData() const;
	Aspect aspect() const;

	void missingNetworkProbability(double probability);
	double missingNetworkProbability() const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability() const;

// The currentPermutationLength defines the length of the stretch
// that is permuted in a **permute step.
// It is between a minimum and a maximum value, defined in the algorithm.
// It is dynamically updated, depending on the acceptance rates.
// That is why currentPermutationLength is not integer.

	void currentPermutationLength(double value);
	double currentPermutationLength() const;

	void updateCurrentPermutationLength(bool accept);

	void createEndStateDifferences();

	void recordOutcome(const MiniStep & miniStep, bool accept,
		int stepType, bool misdat);

    bool smallNeighbourhoodChange(MiniStep * pMiniStep1, MiniStep * pMiniStep2,
                             DependentVariable * pVariable,
                             NetworkVariable * pNetworkVariable,
                             int ego1, int alter1);
    bool neighbourhoodChange(MiniStep * pMiniStep1, MiniStep * pMiniStep2,
                             DependentVariable * pVariable,
                             NetworkVariable * pNetworkVariable,
                             int ego1, int alter1);

	void gotoLastState();


private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();
	bool validInsertMissingStep(const Option * pOption,
		int d0,
		const MiniStep * pMiniStepA);
	bool validDeleteMissingStep(MiniStep * pMiniStepA, bool applyTwice);
	MiniStep * createMiniStep(const Option * pOption,
		int difference, bool value) const;

	double lproposalProbability{};
	bool lmissingData{};
	Aspect laspect{};
	double lprobabilityArray[8]{};// probabilities of MH step types
//	int lacceptances[NBRTYPES];
//	int lrejections[NBRTYPES];
	int laborts[NBRTYPES]{};
	double lmissingNetworkProbability{};
	double lmissingBehaviorProbability{};
	// current length of permuted interval
	double lcurrentPermutationLength{};
	unsigned lthisPermutationLength{};

	// A vector of options with missing values in the initial observation
	std::vector<const Option *> linitialMissingOptions{};
};

}

#endif /* MLSIMULATION_H_ */
