/*
 * SettingInfo.h
 *
 *  Created on: May 6, 2015
 *      Author: ortmann
 */

#ifndef SRC_MODEL_SETTINGS_SETTINGINFO_H_
#define SRC_MODEL_SETTINGS_SETTINGINFO_H_

#include <string>

namespace siena {

class Permission_Type {
public:
	static const Permission_Type BOTH;
	static const Permission_Type UP;
	static const Permission_Type DOWN;

	bool operator==(const Permission_Type& rhs) const {
		return value == rhs.value;
	}

	bool operator!=(const Permission_Type& rhs) const {
		return !operator==(rhs);
	}

	Permission_Type(const Permission_Type &rhs) :
			value(rhs.value) {
	}

	Permission_Type & operator=(const Permission_Type &rhs) {
		value = rhs.value;
		return *this;
	}

private:

	Permission_Type(int val) :
			value(val) {

	}

	int value;
};

class SettingInfo {
public:
	SettingInfo(const std::string& id, const std::string& settingType,
			const std::string& covarName, const Permission_Type permType);
	virtual ~SettingInfo();

	const std::string& getId() const;
	const std::string& getSettingType() const;
	const std::string& getCovarName() const;
	const Permission_Type getPermType() const;

private:

	std::string lid {};

	std::string lsettingType {};

	std::string lcovarName {};

	Permission_Type lpermType;

};

} /* namespace siena */

#endif /* SRC_MODEL_SETTINGS_SETTINGINFO_H_ */
