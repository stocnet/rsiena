/*
 * SettingsFactory.cpp
 *
 *  Created on: 24.06.2014
 *      Author: ortmann
 */

#include <sstream>

#include "SettingsFactory.h"
#include "Setting.h"
#include "UniversalSetting.h"
#include "PrimarySetting.h"
#include "DyadicSetting.h"
#include "MeetingSetting.h"
#include "ComposableSetting.h"
#include <stdexcept>

using namespace std;

namespace siena {

SettingsFactory::SettingsFactory() {
}

SettingsFactory::~SettingsFactory() {
}

Setting* SettingsFactory::createSetting(const SettingInfo& settinginfo) {
	const std::string& settingName = settinginfo.getSettingType();
	if (settingName == "primary") {
		if (settinginfo.getCovarName().empty()) {
		return new PrimarySetting();
		}
		return new ComposableSetting(new PrimarySetting(), new DyadicSetting());
	} else if (settingName == "dyadic") {
		return new MeetingSetting(new DyadicSetting(),
				settinginfo.getPermType());
	} else if (settingName == "universal") {
		return new MeetingSetting(new UniversalSetting(),
				settinginfo.getPermType());
	} else {
		throw invalid_argument("wrong name");
	}
	return 0;
}

std::vector<std::string> SettingsFactory::split(std::string str,
		char delimiter) {
	std::vector<std::string> substrings;
	std::stringstream ss(str);
	std::string tok;
	while (std::getline(ss, tok, delimiter)) {
		substrings.push_back(tok);
	}
	return substrings;
}

} /* namespace siena */
