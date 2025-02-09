/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: https://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LogTable.cpp
 *
 * Description: This file contains the implementation of the LogTable class.
 * It returns log(0) = 0, and further log(i) = natural logarithm of i.
 *****************************************************************************/

#include <cmath>
#include "LogTable.h"

namespace siena
{

// The size of the lookup table. The natural logarithms of integers larger than
// this will be calculated with log() directly.

const int LIMIT = 1000;


// Nullify the pointer to the single instance so we can allocate it
// on demand.

LogTable * LogTable::lpInstance = 0;


/**
 * Returns the only instance of the LogTable class. The instance
 * is allocated on demand.
 */
LogTable * LogTable::instance()
{
	if (!LogTable::lpInstance)
	{
		LogTable::lpInstance = new LogTable();
	}

	return LogTable::lpInstance;
}


/**
 * Constructor.
 */
LogTable::LogTable()
{
	this->ltable = new double[LIMIT];

	// Leave the entries undefined, so they can be computed on demand.

	this->ltable[0] = 0.0;

	for (int i = 1; i < LIMIT; i++)
	{
		this->ltable[i] = -1;
	}
}


/**
 * Destructor.
 */
LogTable::~LogTable()
{
	delete[] this->ltable;
}


/**
 * Returns the natural logarithm of the given integer.
 */
double LogTable::log(int i)
{
	double log = 0;

	if (i < LIMIT)
	{
		if (this->ltable[i] < 0)
		{
			this->ltable[i] = std::log((double) i);
		}

		log = this->ltable[i];
	}
	else
	{
		log = std::log((double) i);
	}

	return log;
}

}
