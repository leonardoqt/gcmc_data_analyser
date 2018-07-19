#include <cstring>
#include "class.h"

atom& atom :: operator=(const atom& B)
{
	symbol = B.symbol;
	sym = B.sym;
	pos = B.pos;
	return *this;
}
