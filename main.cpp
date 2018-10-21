#include <iostream>
#include <cstring>
#include <iomanip>
#include "class.h"

using namespace std;

int main()
{
	ifstream input;
	ofstream output, sdd, image;
	cell c1;
	int num_el=0, num_file=0;
	string *el, file_name;
	double max_bond_length, h_surf;
	double charge[2]={1.0,-2.0};
	
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
	sdd.open("sdd.dat");
	output<<setw(17)<<"file name |"<<setw(9*num_el*num_el)<<"bond length |"<<setw(9*num_el*num_el)<<"coord num |"<<setw(9*num_el)<<"composition |"<<setw(9*num_el*num_el)<<"bond vec norm |"<<setw(9*num_el*num_el)<<"bond vec z |"<<setw(9)<<"gii |"<<setw(9)<<"bv |"<<setw(9)<<"bvv |"<<setw(13)<<"energy"<<endl;
	sdd<<setw(17)<<"file name |"<<setw(9*num_el*num_el)<<"bond length |"<<setw(9*num_el*num_el)<<"coord num |"<<setw(9*num_el*num_el)<<"bond vec norm |"<<setw(9*num_el*num_el)<<"bond vec z"<<endl;

	output<<setw(17)<<' ';
	sdd<<setw(17)<<' ';
	// bond length
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
		{
			output<<setw(9)<<el[t1]+"--"+el[t2];
			sdd<<setw(9)<<el[t1]+"--"+el[t2];
		}
	// coord num
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
		{
			output<<setw(9)<<el[t1]+"--"+el[t2];
			sdd<<setw(9)<<el[t1]+"--"+el[t2];
		}
	// composition
	for(int t1=0; t1<num_el; t1++)
		output<<setw(9)<<el[t1];
	// bond vec norm
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
		{
			output<<setw(9)<<el[t1]+"--"+el[t2];
			sdd<<setw(9)<<el[t1]+"--"+el[t2];
		}
	// bond vec z
	for(int t1=0; t1<num_el; t1++)
		for(int t2=0; t2<num_el; t2++)
		{
			output<<setw(9)<<el[t1]+"--"+el[t2];
			sdd<<setw(9)<<el[t1]+"--"+el[t2];
		}
	output<<endl;
	sdd<<endl;

	for (int t1=0; t1<num_file; t1++)
	{
		cin>>file_name;
		if (cin.eof())
		{
			cout<<"Not enough file !"<<endl;
			exit(0);
		}
		input.open(file_name);
//		image.open(file_name+".csv");
		c1.read_file(input);
		c1.build_e_l(el,num_el);
		c1.find_neighbor(max_bond_length);
		c1.get_mean_bond_length();
		c1.get_mean_coord_num_surf(h_surf);
		c1.get_composition_surf(h_surf);
		c1.get_mean_bond_vec_norm(h_surf);
		// !!!!! not transferrable
		c1.get_gii(h_surf, "Ag", "O",1.842,0.37,1,2);
		c1.get_bv(h_surf, "Ag", "O",1.842,0.37);
		c1.get_bvv(h_surf, "Ag", "O",1.842,0.37);
		output<<setw(17)<<file_name;
		sdd<<setw(17)<<file_name;
		c1.print_all(output);
		c1.print_all_sdd(sdd);
//		c1.generate_image(40,40,image,charge,h_surf);
		c1.clean();
		input.close();
//		image.close();
	}

	output.close();
	sdd.close();
	return 0;
}
