/*
 * SettingInfo.cpp
 *
 *  Created on: May 6, 2015
 *      Author: ortmann
 */

#include <model/settings/SettingInfo.h>

namespace siena {

const Permission_Type Permission_Type::BOTH = Permission_Type(0);
const Permission_Type Permission_Type::UP = Permission_Type(1);
const Permission_Type Permission_Type::DOWN = Permission_Type(2);

SettingInfo::SettingInfo(const std::string& id, const std::string& settingType,
		const std::string& covarName, const Permission_Type permType) :
		lid(id), //
		lsettingType(settingType), //
		lcovarName(covarName), //
		lpermType(permType) {
}

SettingInfo::~SettingInfo() {
}

const std::string& SettingInfo::getId() const {
	return lid;
}

const std::string& SettingInfo::getSettingType() const {
	return lsettingType;
}

const std::string& SettingInfo::getCovarName() const {
	return lcovarName;
}

const Permission_Type SettingInfo::getPermType() const {
	return lpermType;
}

} /* namespace siena */
