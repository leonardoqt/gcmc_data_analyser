#include <iostream>
#include <cstring>
#include "class.h"

using namespace std;

int main()
{
	ifstream input;
	cell c1;
	string el[3]={"Ag","A","g"};

	input.open("structure0001.xsf");
	c1.read_file(input);
	c1.build_e_l(el,3);
	c1.find_neighbor(3.0);
	
	c1.get_mean_bond_length();
	c1.get_mean_coord_num_surf(4.0);
	c1.get_composition_surf(4.0);
	c1.get_mean_bond_vec_norm(4.0);
	c1.clean();

	cout<<endl;

	input.open("structure0001.xsf");
	c1.read_file(input);
	c1.build_e_l(el,3);
	c1.find_neighbor(3.0);
	
	c1.get_mean_bond_length();
	c1.get_mean_coord_num_surf(4.0);
	c1.get_composition_surf(4.0);
	c1.get_mean_bond_vec_norm(4.0);

	return 0;
}
