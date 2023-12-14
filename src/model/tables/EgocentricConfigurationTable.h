#ifndef EGOCENTRICCONFIGURATIONTABLE_H_
#define EGOCENTRICCONFIGURATIONTABLE_H_

#include "ConfigurationTable.h"

namespace siena
{

/**
 * This class defines a table storing a number of network configurations of
 * certain type per actor.
 *
 * There can be various types of configurations (like the number of two-paths,
 * etc.), each implemented in a derived class. Normally, the derived classes
 * should implement only the purely virtual method vCalculate.
 *
 * The configuration tables are grouped in and owned by an instance of the
 * NetworkCache class.
 *
 * The configuration tables are used to speedup the calculations involving
 * effects.
 */
class EgocentricConfigurationTable: public ConfigurationTable
{
public:
	EgocentricConfigurationTable(NetworkCache * pOwner);
	virtual ~EgocentricConfigurationTable();

	void initialize(int ego);
	virtual int get(int i);

protected:
	int ego() const;

private:
	// Indicates the ego this table has been calculated for
	int lego {};

	bool lupdated {};
};

}

#endif /* EGOCENTRICCONFIGURATIONTABLE_H_ */
