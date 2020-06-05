/*
 * SettingsFactory.h
 *
 *  Created on: 24.06.2014
 *      Author: ortmann
 */

#ifndef SETTINGSFACTORY_H_
#define SETTINGSFACTORY_H_

#include <string>
#include <vector>

#include "SettingInfo.h"

namespace siena {

class Setting;

class SettingsFactory {
public:
	SettingsFactory();
	virtual ~SettingsFactory();

	Setting* createSetting(const SettingInfo& info);
private:
	std::vector<std::string> split(std::string str, char delimiter);

};

} /* namespace siena */
#endif /* SETTINGSFACTORY_H_ */
