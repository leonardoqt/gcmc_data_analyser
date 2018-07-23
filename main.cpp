#include <iostream>
#include <cstring>
#include <iomanip>
#include "class.h"

using namespace std;

int main()
{
	ifstream input;
	ofstream output;
	cell c1;
	int num_el=0, num_file=0;
	string *el, file_name;
	double max_bond_length, h_surf;
	
	cin>>max_bond_length>>h_surf;
	cin>>num_el;
	if (num_el == 0)
	{
		cout<<"Number of atom not specified !"<<endl;
		exit(0);
	}

	el = new string[num_el];
	for (int t1=0; t1<num_el; t1++)
		cin>>el[t1];

	cin>>num_file;
	if (num_file == 0)
	{
		cout<<"Number of files not specified !"<<endl;
		exit(0);
	}
	
	output.open("descriptor.dat");
	output<<setw(17)<<"file name |"<<setw(9*num_el*num_el)<<"bond length |"<<setw(9*num_el*num_el)<<"coord num |"<<setw(9*num_el)<<"composition |"<<setw(9*num_el*num_el)<<"bond vec norm |"<<setw(9*num_el*num_el)<<"bond vec z |"<<setw(9)<<"gii |"<<setw(13)<<"energy"<<endl;

	output<<setw(17)<<' ';
	// bond length
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
			output<<setw(9)<<el[t1]+"--"+el[t2];
	// coord num
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
			output<<setw(9)<<el[t1]+"--"+el[t2];
	// composition
	for(int t1=0; t1<num_el; t1++)
		output<<setw(9)<<el[t1];
	// bond vec norm
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
			output<<setw(9)<<el[t1]+"--"+el[t2];
	// bond vec z
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
			output<<setw(9)<<el[t1]+"--"+el[t2];
	output<<endl;

	for (int t1=0; t1<num_file; t1++)
	{
		cin>>file_name;
		if (cin.eof())
		{
			cout<<"Not enough file !"<<endl;
			exit(0);
		}
		input.open(file_name);
		c1.read_file(input);
		c1.build_e_l(el,num_el);
		c1.find_neighbor(max_bond_length);
		c1.get_mean_bond_length();
		c1.get_mean_coord_num_surf(h_surf);
		c1.get_composition_surf(h_surf);
		c1.get_mean_bond_vec_norm(h_surf);
		// !!!!! not transferrable
		c1.get_gii(h_surf, "Ag", "O",1.842,0.37,1,2);
		output<<setw(17)<<file_name;
		c1.print_all(output);
		c1.clean();
		input.close();
	}

/*
	c1.read_file(input);
	c1.build_e_l(el,2);
	c1.find_neighbor(3.0);
	
	c1.get_mean_bond_length();
	c1.get_mean_coord_num_surf(4.0);
	c1.get_composition_surf(4.0);
	c1.get_mean_bond_vec_norm(4.0);
	c1.print_all(out);
	c1.clean();

	input2.open("structure0001.xsf");
	c1.read_file(input2);
	c1.build_e_l(el,2);
	c1.find_neighbor(3.0);
	
	c1.get_mean_bond_length();
	c1.get_mean_coord_num_surf(4.0);
	c1.get_composition_surf(4.0);
	c1.get_mean_bond_vec_norm(4.0);
	c1.print_all(out);
*/
	return 0;
}
