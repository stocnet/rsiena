/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectValueTable.cpp
 *
 * Description: This file contains the implementation of the
 * EffectValueTable class.
 *****************************************************************************/
#include <cmath>

#include "EffectValueTable.h"

namespace siena {

/**
 * Creates a new look-up table for <i>n</i> actors and the given function.
 */
EffectValueTable::EffectValueTable(int n, double (*pFunction)(int)) {
	this->lpFunction = pFunction;
	this->lvalues = new double[n];
	this->lparameterValues = new double[n];

	// Calculate the table for parameter 0

	this->lparameter = 0;

	for (int i = 0; i < n; i++) {
		// exp(0) = 1

		this->lvalues[i] = 1;
		this->lparameterValues[i] = 0;
	}
}

/**
 * Deallocates this look-up table.
 */
EffectValueTable::~EffectValueTable() {
	delete[] this->lvalues;
	delete[] this->lparameterValues;

	this->lvalues = 0;
	this->lparameterValues = 0;
	this->lpFunction = 0;
}

/**
 * Returns the current value of the effect parameter.
 */
double EffectValueTable::parameter() const {
	return this->lparameter;
}

/**
 * Stores the effect parameter.
 */
void EffectValueTable::parameter(double value) {
	this->lparameter = value;
}

/**
 * Returns the value of the effect for the argument <i>i</i>.
 */
double EffectValueTable::value(int i) {
	if (this->lparameterValues[i] != this->lparameter) {
		// The value stored in the table was calculated for a different
		// parameter, hence we must recalculate the value.

		this->lvalues[i] = std::exp(this->lparameter * this->lpFunction(i));
		this->lparameterValues[i] = this->lparameter;
	}

	return this->lvalues[i];
}

}
