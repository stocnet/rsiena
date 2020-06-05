/*
 * PrimarySetting.h
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#ifndef PRIMARYSETTING_H_
#define PRIMARYSETTING_H_

#include "GeneralSetting.h"
#include "network/layers/DistanceTwoLayer.h"
#include "network/layers/PrimaryLayer.h"

namespace siena {

class Network;

class PrimarySetting: public GeneralSetting {
public:
	PrimarySetting();

	virtual ~PrimarySetting();

	/** @copydoc ASetting::initSetting(Network* const lpNetwork) */
	virtual void initSetting(Network* const lpNetwork);
	/** @copydoc ASetting::terminateSetting(Network* const lpNetwork) */
	virtual void terminateSetting(Network* const lpNetwork);

	const Network * pPrimaryNetwork() const;

	ITieIterator* getSteps();

	int getSize();


protected:
	void initSetting(); // called on each initSetting(int ego)
	void terminateSetting();

private:
	Network* lpNetwork; // the underlying dependent

	// old
	// DistanceTwoLayer rDistTwoLayer;
	// ITieIterator* lpiter;

	// new
	PrimaryLayer lPLayer; // the 2-path extension

};

} /* namespace siena */
#endif /* PRIMARYSETTING_H_ */
