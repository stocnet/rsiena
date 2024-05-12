/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateMixedNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateMixedNetworkAlterFunction class.
 *****************************************************************************/

// combines MixedNetworkAlterFunction and CovariateNetworkAlterFunction

#ifndef DYADICCOVARIATEMIXEDNETWORKALTERFUNCTION_H_ 
#define DYADICCOVARIATEMIXEDNETWORKALTERFUNCTION_H_ 

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"
#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class TwoNetworkCache;
class NetworkCache;
class ConstantDyadicCovariate;
class ChangingDyadicCovariate;

class DyadicCovariateMixedNetworkAlterFunction: public MixedNetworkAlterFunction
{
public:
	DyadicCovariateMixedNetworkAlterFunction(std::string firstNetworkName,
		std::string secondNetworkName, std::string dyadicCovariateName);
	virtual ~DyadicCovariateMixedNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

protected:
	double dyadicValue(int i, int j) const;
	bool missing(int i, int j) const;
	bool constantCovariate() const;
	int dyCov_n() const;
	int dyCov_m() const;
	int firstNet_n() const;
	int firstNet_m() const;
	int secondNet_n() const;
	int secondNet_m() const;

private:
	std::string lDyadicCovariateName;
	int lperiod {};
    // number of elements of the respective node sets:
	int ldyadic_n {};
	int ldyadic_m {};
	int lFirstNet_n {};
	int lFirstNet_m {};
	int lSecondNet_n {};
	int lSecondNet_m {};
	// The constant covariate this effect depends on or 0, if the
	// effect depends on a changing covariate:
	ConstantDyadicCovariate * lpConstantCovariate;

	// The changing covariate this effect depends on or 0, if the
	// effect depends on a constant covariate:
	ChangingDyadicCovariate * lpChangingCovariate;

	// flag to control exclusion of missing values:	
	bool lexcludeMissings {};

	/* -Wunused-private-field
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	TwoNetworkCache * lpTwoNetworkCache;
	*/

};


}

#endif /* DYADICCOVARIATEMIXEDNETWORKALTERFUNCTION_H_ */
