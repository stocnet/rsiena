/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConfigurationTable.h
 *
 * Description: This file defines the class ConfigurationTable.
 *****************************************************************************/

#ifndef CONFIGURATIONTABLE_H_
#define CONFIGURATIONTABLE_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkCache;
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
class ConfigurationTable
{
public:
	ConfigurationTable(NetworkCache * pOwner);
	virtual ~ConfigurationTable();

	virtual int get(int i);

protected:
	NetworkCache * pOwner() const;
	const Network * pNetwork() const;

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
	NetworkCache * lpOwner;

	// The network this configuration table is associated with
	const Network * lpNetwork;

	// The modification count of the network on the last time this table
	// was calculated.

	int llastModificationCount;
};

}

#endif /*CONFIGURATIONTABLE_H_*/
