#ifndef SQRTTABLE_H_
#define SQRTTABLE_H_

namespace siena
{

/**
 * A singleton class that provides an efficient computation of square roots
 * of integers. Once a root is calculated, it is stored in a table for
 * subsequent reuse.
 *
 * Usage example:
 * SqrtTable * pTable = SqrtTable::instance();
 * double root2 = pTable->sqrt(2);
 */
class SqrtTable
{
public:
	static SqrtTable * instance();
	double sqrt(int i);

private:
	// The constructor and destructor have to be private such that instances
	// of this class cannot be created directly outside the class' scope.

	SqrtTable();
	virtual ~SqrtTable();

	// The single instance of this class
	static SqrtTable * lpInstance;

	// A table storing the square roots that have been calculated before
	// (ltable[i] = -1 if the root has not been calculated for the number i)

	double * ltable {};
};

}

#endif /*SQRTTABLE_H_*/
