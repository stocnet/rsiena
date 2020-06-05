/*
 * MeetingSetting.h
 *
 *  Created on: 23.06.2014
 *      Author: ortmann
 */

#ifndef MEETINGSETTING_H_
#define MEETINGSETTING_H_

#include "Setting.h"
#include <vector>
#include "SettingInfo.h"

namespace siena {

class MeetingSetting: public Setting {
public:

	virtual ~MeetingSetting();

	MeetingSetting(Setting* setting, const Permission_Type type);

	void initSetting(Network* const lpNetwork);

	void terminateSetting(Network* const lpNetwork);

	void initDyadicSetting(const std::map<int, double>& row, int ego);

	void initPermittedSteps(const bool* const permitted);

	ITieIterator* getSteps();

	int getSize();

	ITieIterator* getPermittedSteps();

	int getPermittedSize();

	virtual bool validate(const Network* const lpNetwork);

protected:

	void initSetting();

	void terminateSetting();

private:

	Setting* lpSetting;

	ITieIterator* lpPermittedSteps;

	const Permission_Type lPermissionType;

};

} /* namespace siena */

#endif /* MEETINGSETTING_H_ */
