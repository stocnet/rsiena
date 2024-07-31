/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * MixedNetworkAlterFunction class.
 *****************************************************************************/


#ifndef MIXEDNETWORKALTERFUNCTION_H_
#define MIXEDNETWORKALTERFUNCTION_H_

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"

namespace siena
{

class Network;
class TwoNetworkCache;
class NetworkCache;
class CommonNeighborIterator;


class MixedNetworkAlterFunction: public AlterFunction
{
public:
	MixedNetworkAlterFunction(std::string firstNetworkName,
		std::string secondNetworkName);
	// for gmom:
	MixedNetworkAlterFunction(std::string firstNetworkName,
		std::string secondNetworkName, const bool simulatedState);
	
	virtual ~MixedNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void initialize(const Data * pData,
		State * pState, State * pSimulatedState, int period, Cache * pCache);

protected:
	bool firstOutTieExists(int alter) const;
	bool secondOutTieExists(int alter) const;
	inline const Network * pFirstNetwork() const;
	inline const Network * pSecondNetwork() const;
	inline TwoNetworkCache * pTwoNetworkCache() const;
	inline NetworkCache * pFirstNetworkCache() const;
	CommonNeighborIterator firstNetworkInStars(int i, int j) const;

private:
	std::string lfirstNetworkName {};
	std::string lsecondNetworkName {};
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	TwoNetworkCache * lpTwoNetworkCache;
	NetworkCache * lpFirstNetworkCache;
	//! If `1` value(), pNetwork() and pNetworkCache() etc. return
	//! the simulated value at the end of the period:
	int lSimulatedOffset {};
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * MixedNetworkAlterFunction::pFirstNetwork() const
{
	return this->lpFirstNetwork;
}

const Network * MixedNetworkAlterFunction::pSecondNetwork() const
{
	return this->lpSecondNetwork;
}

TwoNetworkCache * MixedNetworkAlterFunction::pTwoNetworkCache() const
{
	return this->lpTwoNetworkCache;
}

NetworkCache * MixedNetworkAlterFunction::pFirstNetworkCache() const
{
	return this->lpFirstNetworkCache;
}

}

#endif /* MIXEDNETWORKALTERFUNCTION_H_ */
