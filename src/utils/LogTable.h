#ifndef LOGTABLE_H_
#define LOGTABLE_H_

namespace siena
{

/**
 * A singleton class that provides an efficient computation of square roots
 * of integers. Once a root is calculated, it is stored in a table for
 * subsequent reuse.
 *
 * Usage example:
 * LogTable * pTable = LogTable::instance();
 * double log2 = pTable->log(2);
 */
class LogTable
{
public:
	static LogTable * instance();
	double log(int i);

private:
	// The constructor and destructor have to be private such that instances
	// of this class cannot be created directly outside the class' scope.

	LogTable();
	virtual ~LogTable();

	// The single instance of this class
	static LogTable * lpInstance;
	// A table storing the logs that have been calculated before;
	// (ltable[i] = -1 if the log has not been calculated for the number i)

	double * ltable {};
};

}

#endif /*LOGTABLE_H_*/
