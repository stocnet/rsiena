/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedConfigurationTable.h
 *
 * Description: This file defines the class MixedConfigurationTable.
 *****************************************************************************/

#ifndef MIXEDCONFIGURATIONTABLE_H_
#define MIXEDCONFIGURATIONTABLE_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class TwoNetworkCache;
class Network;


// ----------------------------------------------------------------------------
// Section: Class description
// ----------------------------------------------------------------------------

/**
 * This class defines a table storing a number of network configurations of
 * certain type per actor.
 *
 * There can be various types of configurations, each implemented in a derived
 * class. Normally, the derived classes should implement only the purely
 * virtual method vCalculate.
 *
 * The configuration tables are grouped in and owned by an instance of the
 * NetworkCache class.
 *
 * The configuration tables are used to speedup the calculations involving
 * effects.
 */
class MixedConfigurationTable
{
public:
	MixedConfigurationTable(TwoNetworkCache * pOwner);
	virtual ~MixedConfigurationTable();

	virtual int get(int i);

protected:
	TwoNetworkCache * pOwner() const;
	const Network * pFirstNetwork() const;
	const Network * pSecondNetwork() const;

	/**
	 * An abstract method that has to be implemented by derived classes
	 * to actually calculate the number of configurations of the respective
	 * type.
	 */
	virtual void calculate() = 0;

	void reset();

	// The internal storage
	int * ltable;

private:
	// The network cache owning this configuration table
	TwoNetworkCache * lpOwner;

	// The first network this configuration table is associated with
	const Network * lpFirstNetwork;

	// The secondnetwork this configuration table is associated with
	const Network * lpSecondNetwork;

	// The modification count of the first network on the last time this table
	// was calculated.

	int llastFirstModificationCount {};

	// The modification count of the second network on the last time this table
	// was calculated.

	int llastSecondModificationCount {};

	// The number of elements in ltable
	int ltableSize {};
};

}

#endif /*MIXEDCONFIGURATIONTABLE_H_*/
