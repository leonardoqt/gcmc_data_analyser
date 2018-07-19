#include <iostream>
#include "class.h"

using namespace std;

int main()
{
	ifstream input;
	cell c1;
	input.open("structure0001.xsf");
	c1.read_file(input);
	c1.find_neighbor(3.0);
	return 0;
}
