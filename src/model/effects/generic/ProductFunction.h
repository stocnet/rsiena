#ifndef PRODUCTFUNCTION_H_
#define PRODUCTFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class ProductFunction: public AlterFunction
{
public:
	ProductFunction(AlterFunction * pFirstFunction,
		AlterFunction * pSecondFunction);
	virtual ~ProductFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double value(int alter);

private:
	AlterFunction * lpFirstFunction;
	AlterFunction * lpSecondFunction;
};

}

#endif /* PRODUCTFUNCTION_H_ */
