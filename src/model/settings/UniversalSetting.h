/*
 * UniversalSetting.h
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#ifndef UNIVERSALSETTING_H_
#define UNIVERSALSETTING_H_

#include "GeneralSetting.h"
#include <vector>

namespace siena {

class UniversalSetting: public GeneralSetting {
public:
	UniversalSetting();

	virtual ~UniversalSetting();

	/**
	 * @copydoc ASetting::initSetting(Network* const lpNetwork)
	 */
	virtual void initSetting(Network* const lpNetwork);

	/**
	 * @copydoc ASetting::terminateSetting(Network* const lpNetwork)
	 */
	virtual void terminateSetting(Network* const lpNetwork);
	using GeneralSetting::terminateSetting;

	int getSize();

	ITieIterator* getSteps();

protected:

	void initSetting();

private:
	std::vector<int> rSteps;

};

} /* namespace siena */
#endif /* UNIVERSALSETTING_H_ */
