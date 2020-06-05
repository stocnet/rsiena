/*
 * DyadicSetting.h
 *
 *  Created on: 24.07.2014
 *      Author: ortmann
 */

#ifndef DYADICSETTING_H_
#define DYADICSETTING_H_

#include "GeneralSetting.h"

namespace siena {

class DyadicSetting: public GeneralSetting {

public:

	DyadicSetting();

	virtual ~DyadicSetting();

	/**
	 * @copydoc ASetting::initSetting(Network* const lpNetwork)
	 */
	void initSetting(Network* const lpNetwork);

	/**
	 * @copydoc ASetting::terminateSetting(Network* const lpNetwork)
	 */
	void terminateSetting(Network* const lpNetwork);

	int getSize();

	ITieIterator* getSteps();

	void initDyadicSetting(const std::map<int, double>& row, int ego);

protected:

	void initSetting();

	void terminateSetting();

private:

	void createSteps(const std::map<int, double>& row);

	ITieIterator* lpiter;

};

} /* namespace siena */
#endif /* DYADICSETTING_H_ */
