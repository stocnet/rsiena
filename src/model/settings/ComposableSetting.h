/*
 * ComposableSetting.h
 *
 *  Created on: May 5, 2015
 *      Author: ortmann
 */

#ifndef SRC_MODEL_SETTINGS_COMPOSABLESETTING_H_
#define SRC_MODEL_SETTINGS_COMPOSABLESETTING_H_

#include "Setting.h"

namespace siena {

class ComposableSetting: public Setting {
public:
	virtual ~ComposableSetting();

	ComposableSetting(Setting* firstSetting, Setting* secondSetting);

	void initSetting(Network* const lpNetwork);

	void terminateSetting(Network* const lpNetwork);

	void initDyadicSetting(const std::map<int, double>& row, int ego);

	void initPermittedSteps(const bool* const permitted);

	ITieIterator* getSteps();

	int getSize();

	ITieIterator* getPermittedSteps();

	int getPermittedSize();

protected:

	void initSetting();

	void terminateSetting();

private:

	Setting* lpfirstSetting;

	Setting* lpsecondSetting;

	ITieIterator* lpSteps;

	ITieIterator* lpPermittedSteps;

};

} /* namespace siena */

#endif /* SRC_MODEL_SETTINGS_COMPOSABLESETTING_H_ */
