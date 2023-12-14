/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedEgocentricConfigurationTable.h
 *
 * Description: This file defines the class MixedEgocentriConfigurationTable.
 *****************************************************************************/
#ifndef MIXEDEGOCENTRICCONFIGURATIONTABLE_H_
#define MIXEDEGOCENTRICCONFIGURATIONTABLE_H_

#include "MixedConfigurationTable.h"

namespace siena
{

/**
 * This class defines a table storing a number of network configurations of
 * certain type per actor.
 *
 * There can be various types of configurations (like the number of two-paths,
 * etc.), each implemented in a derived class. Normally, the derived classes
 * should implement only the purely virtual method vCalculate.
 *
 * The configuration tables are grouped in and owned by an instance of the
 * NetworkCache class.
 *
 * The configuration tables are used to speedup the calculations involving
 * effects.
 */
class MixedEgocentricConfigurationTable: public MixedConfigurationTable
{
public:
	MixedEgocentricConfigurationTable(TwoNetworkCache * pOwner);
	virtual ~MixedEgocentricConfigurationTable();

	void initialize(int ego);
	virtual int get(int i);

protected:
	int ego() const;

private:
	// Indicates the ego this table has been calculated for
	int lego {};

	bool lupdated {};
};

}

#endif /* MIXEDEGOCENTRICCONFIGURATIONTABLE_H_ */
