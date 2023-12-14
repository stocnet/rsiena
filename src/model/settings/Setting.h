/*
 * Setting.h
 *
 *  Created on: 23.07.2014
 *      Author: ortmann
 */

#ifndef SETTING_H_
#define SETTING_H_

#include <map>

namespace siena {

class ITieIterator;
class Network;

class Setting {
public:
	virtual ~Setting() {
	}

	/**
	 * Sets the rate of the setting.
	 * @param[in] rate the new rate
	 */
	void setRate(double rate) {
		lrate = rate;
	}

	/**
	 * Returns the rate of the setting.
	 */
	double getRate() {
		return lrate;
	}

	/**
	 * Adds the view to the network.
	 * @param[in] lpNetwork the network
	 */
	virtual void initSetting(Network* const lpNetwork) =0;

	/**
	 * Removes the view to the network.
	 * @param[in] lpNetwork the network
	 */
	virtual void terminateSetting(Network* const lpNetwork) =0;

	void initSetting(int ego);

	void terminateSetting(int ego);

	virtual void initDyadicSetting(const std::map<int, double>& /*row*/, int /*ego*/) {
	}

	virtual void initPermittedSteps(const bool* const permitted) = 0;

	/**
	 * Returns the candidate set of alters for the given actor.
	 *  NOTE: Delete iterator
	 * @param[in] actor the actor whose candidate set has to be returned
	 */
	virtual ITieIterator* getSteps() = 0;

	virtual int getSize() = 0;

	virtual ITieIterator* getPermittedSteps() = 0;

	virtual int getPermittedSize() = 0;

	virtual bool validate(const Network* const lpNetwork);

protected:

	/* C'stor **/
	Setting() :
			lrate(0), //
			lego(-1) {
	}

	virtual void initSetting() = 0;

	virtual void terminateSetting() = 0;

	inline int ego() {
		return lego;
	}

private:
	/* settings rate **/
	double lrate {};

	int lego {};

};
} /* namespace siena */
#endif /* SETTING_H_ */
