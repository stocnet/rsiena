/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessTable.h
 *
 * Description: This file defines the class BetweennessTable.
 *****************************************************************************/

#ifndef BETWEENNESSTABLE_H_
#define BETWEENNESSTABLE_H_

#include "ConfigurationTable.h"
#include "network/NetworkUtils.h"

namespace siena
{

/**
 * This class defines a table of betweenness values. The betweenness of an
 * actor i is the number of ordered actor pairs (j,h) such that there are
 * ties (j,i) and (i,h), but no tie (j,h).
 */
class BetweennessTable : public ConfigurationTable
{
public:
	BetweennessTable(NetworkCache * pOwner);

protected:
	virtual void calculate();
};

}

#endif /*BETWEENNESSTABLE_H_*/
