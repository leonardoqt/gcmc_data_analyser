#include <iostream>
#include <cstring>
#include <fstream>
#include "class.h"

using namespace std;

void cell :: read_file(ifstream &input)
{
	string tmp;
	string unique_el;
	string label_vec="PRIMVEC", label_coord="PRIMCOORD";
	double v1[3];
	int t1,t2,current;
	// rewind file pointer
	input.clear();
	input.seekg(0);
	// find total energy
	input>>tmp>>tmp>>tmp>>tmp>>tot_energy;
	// find lattice const.
	getline(input,tmp);
	while(tmp.find(label_vec) == string::npos)
		getline(input,tmp);
	input>>v1[0]>>v1[1]>>v1[2];
	getline(input,tmp);
	lattice[0] = v1;
	input>>v1[0]>>v1[1]>>v1[2];
	getline(input,tmp);
	lattice[1] = v1;
	input>>v1[0]>>v1[1]>>v1[2];
	getline(input,tmp);
	lattice[2] = v1;
	// find number of atoms
	while(tmp.find(label_coord) == string::npos)
		getline(input,tmp);
	input>>num_atom;
	getline(input, tmp);
	// define number of atom related array
	a_l = new atom[num_atom];
	neighbor_l = new atom*[num_atom];
	num_neighbor = new int[num_atom];
	for (t1=0; t1<num_atom; t1++)
		neighbor_l[t1] = new atom[max_num_neighbor];
	// read each atom and find unique element
	for (t1=0,unique_el=' ',num_element=0; t1<num_atom; t1++)
	{
		input>>a_l[t1].symbol>>v1[0]>>v1[1]>>v1[2];
		getline(input,tmp);
		a_l[t1].pos = v1;
		if (unique_el.find(' '+a_l[t1].symbol+' ') == string::npos)
		{
			num_element++;
			unique_el = unique_el + ' '+a_l[t1].symbol+' ';
		}
	}
	// define number of element related array
	e_l = new string[num_element];
	neighbor_dis = new double*[num_element];
	mean_bond_length = new double*[num_element];
	mean_coord_num_surf = new double*[num_element];
	for(t1=0;t1<num_element;t1++)
	{
		neighbor_dis[t1] = new double[num_element];
		mean_bond_length[t1] = new double[num_element];
		mean_coord_num_surf[t1] = new double[num_element];
	}
	// find symbol of each element and assign element number to each atom
	// !!!!! Should define absolute order for each element in order to 
	//       make sure everything is correct when rearrangeing atoms in xsf.
	for(t1=current=0;t1<num_atom;t1++)
	{
		for(t2=0;t2<current;t2++)
		{
			if (e_l[t2] == a_l[t1].symbol)
				break;
		}
		if (t2 < current)
			a_l[t1].sym = t2;
		else
		{
			e_l[current] = a_l[t1].symbol;
			a_l[t1].sym = current;
			current++;
		}
	}
/*
	for(t1=0;t1<num_element;t1++)
		cout<<e_l[t1]<<'\t'<<t1<<endl;
	for(t1=0;t1<num_atom;t1++)
		cout<<a_l[t1].symbol<<'\t'<<a_l[t1].sym<<endl;
*/
}

void cell :: find_neighbor(double dis)
{
	// get distance table
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			neighbor_dis[t1][t2] = dis;
	// find neighbors
	for (int t1=0; t1<num_atom; t1++)
	{
		num_neighbor[t1] = 0;
		for(int t2=0; t2<num_atom; t2++)
		{
			if (t2 != t1)
				for(int dx=-1;dx<=1;dx++)
				for(int dy=-1;dy<=1;dy++)
				for(int dz=-1;dz<=1;dz++)
				{
					if (((a_l[t2].pos+lattice[0]*dx+lattice[1]*dy+lattice[2]*dz)-a_l[t1].pos).norm() < neighbor_dis[a_l[t1].sym][a_l[t2].sym])
					{
						neighbor_l[t1][num_neighbor[t1]] = a_l[t2];
						neighbor_l[t1][num_neighbor[t1]].pos = a_l[t2].pos+lattice[0]*dx+lattice[1]*dy+lattice[2]*dz;
						num_neighbor[t1]++;
					}
				}
		}
	}
/*
	for (int t1=0; t1<num_atom; t1++)
	{
		cout<<num_neighbor[t1]<<"  "<<1<<endl;
		for(int t2=0; t2<num_neighbor[t1];t2++)
		{
			cout<<neighbor_l[t1][t2].symbol<<'\t'<<neighbor_l[t1][t2].pos.x[0]<<'\t'<<neighbor_l[t1][t2].pos.x[1]<<'\t'<<neighbor_l[t1][t2].pos.x[2]<<endl;
		}
	}
*/
}

void cell :: find_neighbor(double **dis)
{
	// get distance table
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			neighbor_dis[t1][t2] = dis[t1][t2];
	// find neighbors
	for (int t1=0; t1<num_atom; t1++)
	{
		num_neighbor[t1] = 0;
		for(int t2=0; t2<num_atom; t2++)
		{
			if (t2 != t1)
				for(int dx=-1;dx<=1;dx++)
				for(int dy=-1;dy<=1;dy++)
				for(int dz=-1;dz<=1;dz++)
				{
					if (((a_l[t2].pos+lattice[0]*dx+lattice[1]*dy+lattice[2]*dz)-a_l[t1].pos).norm() < neighbor_dis[a_l[t1].sym][a_l[t2].sym])
					{
						neighbor_l[t1][num_neighbor[t1]] = a_l[t2];
						neighbor_l[t1][num_neighbor[t1]].pos = a_l[t2].pos+lattice[0]*dx+lattice[1]*dy+lattice[2]*dz;
						num_neighbor[t1]++;
					}
				}
		}
	}
/*
	for (int t1=0; t1<num_atom; t1++)
	{
		cout<<num_neighbor[t1]<<"  "<<1<<endl;
		for(int t2=0; t2<num_neighbor[t1];t2++)
		{
			cout<<neighbor_l[t1][t2].symbol<<'\t'<<neighbor_l[t1][t2].pos.x[0]<<'\t'<<neighbor_l[t1][t2].pos.x[1]<<'\t'<<neighbor_l[t1][t2].pos.x[2]<<endl;
		}
	}
*/
}

void cell :: get_mean_bond_length()
{
	int **num_bonds;
	num_bonds = new int*[num_element];
	for(int t1=0; t1<num_element; t1++)
		num_bonds[t1] = new int[num_element];
	
	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
		{
			mean_bond_length[t1][t2] = 0;
			num_bonds[t1][t2] = 0;
		}

	for(int t1=0; t1<num_atom; t1++)
		for(int t2=0; t2<num_neighbor[t1]; t2++)
		{
			num_bonds[a_l[t1].sym][neighbor_l[t1][t2].sym]++;
			mean_bond_length[a_l[t1].sym][neighbor_l[t1][t2].sym] += (neighbor_l[t1][t2].pos-a_l[t1].pos).norm();
		}
	
	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
		{
			if (num_bonds[t1][t2] > 0)
			{
				mean_bond_length[t1][t2] /= num_bonds[t1][t2];
//				cout<<e_l[t1]<<"--"<<e_l[t2]<<": "<<mean_bond_length[t1][t2]<<endl;
			}
		}
}

void cell :: get_mean_coord_num_surf(double h_surf)
{
	int *num_atoms;
	num_atoms = new int[num_element];
	for(int t1=0; t1<num_element; t1++)
	{
		num_atoms[t1]=0;
		for(int t2=0; t2<num_element; t2++)
			mean_coord_num_surf[t1][t2];
	}

	for(int t1=0; t1<num_atom; t1++)
	{
		if (a_l[t1].pos.x[2] > h_surf)
		{
			num_atoms[a_l[t1].sym]++;
			for(int t2=0; t2<num_neighbor[t1]; t2++)
				mean_coord_num_surf[a_l[t1].sym][neighbor_l[t1][t2].sym]++;
		}
	}

	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
			if (num_atoms[t1] > 0)
			{
				mean_coord_num_surf[t1][t2] /= num_atoms[t1];
//				cout<<e_l[t1]<<"--"<<e_l[t2]<<": "<<mean_coord_num_surf[t1][t2]<<endl;
			}
}
