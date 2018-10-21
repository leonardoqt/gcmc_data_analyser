#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "class.h"

using namespace std;

void cell :: clean()
{
	// num_element related variable
	for (int t1=0; t1<num_element; t1++)
	{
		delete[] mean_bond_length[t1];
		delete[] mean_coord_num_surf[t1];
		delete[] mean_bond_vec_norm[t1];
		delete[] mean_bond_vec_z[t1];
		delete[] neighbor_dis[t1];
		// sdd
		delete[] mean_bond_length_sdd[t1];
		delete[] mean_coord_num_surf_sdd[t1];
		delete[] mean_bond_vec_norm_sdd[t1];
		delete[] mean_bond_vec_z_sdd[t1];
	}
	// num_atom related variable
	for (int t1=0; t1<num_atom; t1++)
	{
		delete[] neighbor_l[t1];
	}
	delete[] a_l;
	delete[] e_l;
	delete[] num_neighbor;
	delete[] neighbor_l;
	delete[] neighbor_dis;
	delete[] mean_bond_length;
	delete[] mean_coord_num_surf;
	delete[] composition_surf;
	delete[] mean_bond_vec_norm;
	delete[] mean_bond_vec_z;
	//sdd
	delete[] mean_bond_length_sdd;
	delete[] mean_coord_num_surf_sdd;
	delete[] mean_bond_vec_norm_sdd;
	delete[] mean_bond_vec_z_sdd;
}

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
}

void cell :: build_e_l(string * all_el, int num_el)
{
	int t1,t2;
	num_element = num_el;
	// define number of element related array
	e_l = new string[num_element];
	neighbor_dis = new double*[num_element];
	mean_bond_length = new double*[num_element];
	mean_coord_num_surf = new double*[num_element];
	composition_surf = new double[num_element];
	mean_bond_vec_norm = new double*[num_element];
	mean_bond_vec_z = new double*[num_element];
	//sdd
	mean_bond_length_sdd = new double*[num_element];
	mean_coord_num_surf_sdd = new double*[num_element];
	mean_bond_vec_norm_sdd = new double*[num_element];
	mean_bond_vec_z_sdd = new double*[num_element];
	for(t1=0;t1<num_element;t1++)
	{
		neighbor_dis[t1] = new double[num_element];
		mean_bond_length[t1] = new double[num_element];
		mean_coord_num_surf[t1] = new double[num_element];
		mean_bond_vec_norm[t1] = new double[num_element];
		mean_bond_vec_z[t1] = new double[num_element];
		//sdd
		mean_bond_length_sdd[t1] = new double[num_element];
		mean_coord_num_surf_sdd[t1] = new double[num_element];
		mean_bond_vec_norm_sdd[t1] = new double[num_element];
		mean_bond_vec_z_sdd[t1] = new double[num_element];
	}

	for(t1=0; t1<num_element; t1++)
		e_l[t1] = all_el[t1];
	// find element number to each atom
	for(t1=0;t1<num_atom;t1++)
	{
		for(t2=0;t2<num_element;t2++)
		{
			if (e_l[t2] == a_l[t1].symbol)
				break;
		}
		if (t2 < num_element)
		{
			a_l[t1].sym = t2;
		}
		else
		{
			cout<<"element "<<a_l[t1].symbol<<" not found"<<endl;
			exit(0);
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
			mean_bond_length_sdd[t1][t2] = 0;
			num_bonds[t1][t2] = 0;
		}

	for(int t1=0; t1<num_atom; t1++)
		for(int t2=0; t2<num_neighbor[t1]; t2++)
		{
			num_bonds[a_l[t1].sym][neighbor_l[t1][t2].sym]++;
			mean_bond_length[a_l[t1].sym][neighbor_l[t1][t2].sym] += (neighbor_l[t1][t2].pos-a_l[t1].pos).norm();
			mean_bond_length_sdd[a_l[t1].sym][neighbor_l[t1][t2].sym] += pow((neighbor_l[t1][t2].pos-a_l[t1].pos).norm(),2);
		}
	
	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
		{
			if (num_bonds[t1][t2] > 0)
			{
				mean_bond_length[t1][t2] /= num_bonds[t1][t2];
				mean_bond_length_sdd[t1][t2] /= num_bonds[t1][t2];
				mean_bond_length_sdd[t1][t2] = sqrt(mean_bond_length_sdd[t1][t2]-pow(mean_bond_length[t1][t2],2));
//				cout<<e_l[t1]<<"--"<<e_l[t2]<<": "<<mean_bond_length[t1][t2]<<endl;
			}
		}

	for(int t1=0; t1<num_element; t1++)
		delete[] num_bonds[t1];
	delete[] num_bonds;
}

void cell :: get_mean_coord_num_surf(double h_surf)
{
	int *num_atoms, *num_atoms_nei;
	num_atoms = new int[num_element];
	num_atoms_nei = new int[num_element];
	for(int t1=0; t1<num_element; t1++)
	{
		num_atoms[t1]=0;
		for(int t2=0; t2<num_element; t2++)
		{
			mean_coord_num_surf[t1][t2]=0;
			mean_coord_num_surf_sdd[t1][t2]=0;
		}
	}

	for(int t1=0; t1<num_atom; t1++)
	{
		if (a_l[t1].pos.x[2] > h_surf)
		{
			// initiate num_atoms_nei
			for (int t2=0; t2<num_element; t2++)
				num_atoms_nei[t2] = 0;
			// add number of centra atom
			num_atoms[a_l[t1].sym]++;
			// count number of neighbors
			for (int t2=0; t2<num_neighbor[t1]; t2++)
				num_atoms_nei[neighbor_l[t1][t2].sym]++;
			// add to mean and sdd
			for (int t2=0; t2<num_element; t2++)
			{
				mean_coord_num_surf[a_l[t1].sym][t2] += num_atoms_nei[t2];
				mean_coord_num_surf_sdd[a_l[t1].sym][t2] += pow(num_atoms_nei[t2],2);
			}
		}
	}

	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
			if (num_atoms[t1] > 0)
			{
				mean_coord_num_surf[t1][t2] /= num_atoms[t1];
				mean_coord_num_surf_sdd[t1][t2] /= num_atoms[t1];
				mean_coord_num_surf_sdd[t1][t2] = sqrt(mean_coord_num_surf_sdd[t1][t2] - pow(mean_coord_num_surf[t1][t2],2));
//				cout<<e_l[t1]<<"--"<<e_l[t2]<<": "<<mean_coord_num_surf[t1][t2]<<endl;
			}

	delete[] num_atoms;
	delete[] num_atoms_nei;
}

void cell :: get_composition_surf(double h_surf)
{
	double tot;

	for (int t1=0; t1<num_element; t1++)
		composition_surf[t1] = 0;

	for (int t1=0; t1<num_atom; t1++)
		if (a_l[t1].pos.x[2] > h_surf)
			composition_surf[a_l[t1].sym]++;

	tot = 0;
	for (int t1=0; t1<num_element; t1++)
		tot += composition_surf[t1];
	
	for (int t1=0; t1<num_element; t1++)
	{
		composition_surf[t1] /= tot;
//		cout<<e_l[t1]<<": "<<composition_surf[t1]<<endl;
	}
}

void cell :: get_mean_bond_vec_norm(double h_surf)
{
	int **num_bond_vec;
	vec *tmp;
	num_bond_vec = new int*[num_element];
	tmp = new vec[num_element];
	for(int t1=0; t1<num_element; t1++)
	{
		num_bond_vec[t1] = new int[num_element];
		tmp[t1].clean();
	}
	
	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
		{
			mean_bond_vec_norm[t1][t2] = 0;
			mean_bond_vec_z[t1][t2] = 0;
			num_bond_vec[t1][t2] = 0;
			//sdd
			mean_bond_vec_norm_sdd[t1][t2] = 0;
			mean_bond_vec_z_sdd[t1][t2] = 0;
		}

	for(int t1=0; t1<num_atom; t1++)
		if (a_l[t1].pos.x[2] > h_surf)
		{
			for(int t2=0; t2<num_element; t2++)
			{
				tmp[t2].clean();
			}
			for(int t2=0; t2<num_neighbor[t1]; t2++)
			{
				tmp[neighbor_l[t1][t2].sym] = tmp[neighbor_l[t1][t2].sym] + (neighbor_l[t1][t2].pos-a_l[t1].pos);
			}
			for(int t2=0; t2<num_element; t2++)
			{
				mean_bond_vec_norm[a_l[t1].sym][t2] += tmp[t2].norm();
				mean_bond_vec_z[a_l[t1].sym][t2] += tmp[t2].x[2];
				mean_bond_vec_norm_sdd[a_l[t1].sym][t2] += pow(tmp[t2].norm(),2);
				mean_bond_vec_z_sdd[a_l[t1].sym][t2] += pow(tmp[t2].x[2],2);
				if (tmp[t2].norm() > 0)
					num_bond_vec[a_l[t1].sym][t2]++;
			}
		}
	
	for(int t1=0; t1<num_element; t1++)
		for(int t2=0; t2<num_element; t2++)
		{
			if (num_bond_vec[t1][t2] > 0)
			{
				mean_bond_vec_norm[t1][t2] /= num_bond_vec[t1][t2];
				mean_bond_vec_z[t1][t2] /= num_bond_vec[t1][t2];
				mean_bond_vec_norm_sdd[t1][t2] /= num_bond_vec[t1][t2];
				mean_bond_vec_norm_sdd[t1][t2] = sqrt(mean_bond_vec_norm_sdd[t1][t2] - pow(mean_bond_vec_norm[t1][t2],2));
				mean_bond_vec_z_sdd[t1][t2] /= num_bond_vec[t1][t2];
				mean_bond_vec_z_sdd[t1][t2] = sqrt(mean_bond_vec_z_sdd[t1][t2] - pow(mean_bond_vec_z[t1][t2],2));
//				cout<<e_l[t1]<<"--"<<e_l[t2]<<": "<<mean_bond_vec_norm[t1][t2]<<endl;
			}
		}

	for(int t1=0; t1<num_element; t1++)
		delete[] num_bond_vec[t1];
	delete[] tmp;
	delete[] num_bond_vec;
}

// !!!!! this function is not transferrable
void cell :: get_gii(double h_surf, string el1, string el2, double r0, double cc, double n_ox1, double n_ox2)
{
	int i_el1=-1, i_el2=-1;
	int tmp_num, tmp_num_nei;
	double tmp_vec;
	int t1,t2;
	for (t1=0; t1<num_element; t1++)
	{
		if (e_l[t1] == el1)
			i_el1 = t1;
		if (e_l[t1] == el2)
			i_el2 = t1;
	}
	if (i_el1==-1 or i_el2 ==-1)
	{
		cout<<"Error in calculating global instability index, no such element!"<<endl;
		exit(0);
	}

	gii = 0;
	tmp_num = 0;
	// loop for element 1
	for (t1=0; t1<num_atom; t1++)
		if (a_l[t1].sym == i_el1 and a_l[t1].pos.x[2] >= h_surf)
		{
			tmp_num_nei = 0;
			tmp_vec = 0;
			for (t2=0; t2<num_neighbor[t1]; t2++)
				if (neighbor_l[t1][t2].sym == i_el2)
				{
					tmp_vec += pow( (neighbor_l[t1][t2].pos-a_l[t1].pos).norm()/r0,cc );
					tmp_num_nei++;
				}
			if (tmp_num_nei >0)
			{
				gii += (tmp_vec - n_ox1) * (tmp_vec - n_ox1);
				tmp_num++;
			}
		}

	// loop for element 2
	for (t1=0; t1<num_atom; t1++)
		if (a_l[t1].sym == i_el2 and a_l[t1].pos.x[2] >= h_surf)
		{
			tmp_num_nei = 0;
			tmp_vec = 0;
			for (t2=0; t2<num_neighbor[t1]; t2++)
				if (neighbor_l[t1][t2].sym == i_el1)
				{
					tmp_vec += pow( (neighbor_l[t1][t2].pos-a_l[t1].pos).norm()/r0,cc );
					tmp_num_nei++;
				}
			if (tmp_num_nei >0)
			{
				gii += (tmp_vec - n_ox2) * (tmp_vec - n_ox2);
				tmp_num++;
			}
		}
	if (tmp_num >0)
		gii = sqrt(gii / tmp_num);
}

// !!!!! this function is not transferrable
void cell :: get_bv(double h_surf, string el1, string el2, double r0, double cc)
{
	int i_el1=-1, i_el2=-1;
	int tmp_num, tmp_num_nei;
	double tmp_vec;
	int t1,t2;
	for (t1=0; t1<num_element; t1++)
	{
		if (e_l[t1] == el1)
			i_el1 = t1;
		if (e_l[t1] == el2)
			i_el2 = t1;
	}
	if (i_el1==-1 or i_el2 ==-1)
	{
		cout<<"Error in calculating bv, no such element!"<<endl;
		exit(0);
	}

	bv = 0;
	tmp_num = 0;
	// loop for element 1
	for (t1=0; t1<num_atom; t1++)
		if (a_l[t1].sym == i_el1 and a_l[t1].pos.x[2] >= h_surf)
		{
			tmp_num++;
			for (t2=0; t2<num_atom; t2++)
				if (a_l[t2].sym == i_el2)
				{
					bv += exp((r0-(a_l[t2].pos-a_l[t1].pos).norm())/cc);
//					bv += pow(r0 / (a_l[t2].pos-a_l[t1].pos).norm(), cc);
				}
		}
	if (tmp_num > 0)
		bv /= tmp_num;
}

// !!!!! this function is not transferrable
void cell :: get_bvv(double h_surf, string el1, string el2, double r0, double cc)
{
	int i_el1=-1, i_el2=-1;
	int tmp_num, tmp_num_nei;
	vec tmp_vec;
	int t1,t2;
	for (t1=0; t1<num_element; t1++)
	{
		if (e_l[t1] == el1)
			i_el1 = t1;
		if (e_l[t1] == el2)
			i_el2 = t1;
	}
	if (i_el1==-1 or i_el2 ==-1)
	{
		cout<<"Error in calculating bv, no such element!"<<endl;
		exit(0);
	}

	bvv = 0;
	tmp_num = 0;
	// loop for element 1
	for (t1=0; t1<num_atom; t1++)
		if (a_l[t1].sym == i_el1 and a_l[t1].pos.x[2] >= h_surf)
		{
			tmp_num++;
			tmp_vec.clean();
			for (t2=0; t2<num_atom; t2++)
				if (a_l[t2].sym == i_el2)
				{
//					 tmp_vec = tmp_vec + (a_l[t1].pos-a_l[t2].pos)*pow(r0 / (a_l[t2].pos-a_l[t1].pos).norm(), cc);
					 tmp_vec = tmp_vec + (a_l[t1].pos-a_l[t2].pos)*exp((r0-(a_l[t2].pos-a_l[t1].pos).norm())/cc);
				}
			bvv += tmp_vec.norm();
		}
	if (tmp_num > 0)
		bvv /= tmp_num;
}

// not transferrable in terms of defining chagre
void cell :: generate_image(int nx, int ny, std::ofstream& output, double *charge, double h_surf)
{
	double *dx, *dy, **image;
	double height;
	int num_surf;
	double scale = 1.0;
	vec tmp;
	
	dx = new double[nx];
	dy = new double[ny];
	image = new double*[nx];
	for(int t1=0; t1<nx; t1++)
		image[t1] = new double[ny];

	for(int t1=0; t1<nx; t1++)
		for(int t2=0; t2<ny; t2++)
			image[t1][t2] = 0;

	//define dx, dy array
	for(int t1=0; t1<nx; t1++)
		dx[t1] = t1/(double)nx;
	for(int t1=0; t1<ny; t1++)
		dy[t1] = t1/(double)ny;

	// find height as average height
	num_surf = 0;
	height = 0;
	for(int t1=0; t1<num_atom; t1++)
		if (a_l[t1].pos.x[2] > h_surf)
		{
			height += a_l[t1].pos.x[2];
			num_surf++;
		}
	height /= num_surf;
	
	// get image
	for(int t1=0; t1<nx; t1++)
		for(int t2=0; t2<ny; t2++)
		{
			tmp = lattice[0]*dx[t1] + lattice[1]*dy[t2];
			tmp.x[2] = height;
			for(int t3=0; t3<num_atom; t3++)
				if(a_l[t3].pos.x[2]>h_surf)
				{
					for(int aaa=-1; aaa<2; aaa++)
						for(int bbb=-1; bbb<2; bbb++)
						{
							image[t1][t2] += charge[a_l[t3].sym]*exp(-fabs(charge[a_l[t3].sym]) * (tmp - a_l[t3].pos - lattice[0]*aaa - lattice[1]*bbb).norm());
						}
				}
		}
	
	// print out image
	for(int t1=0; t1<nx; t1++)
	{
		output<<setprecision(4)<<fixed<<image[t1][0];
		for(int t2=1; t2<ny; t2++)
			output<<','<<setprecision(4)<<fixed<<image[t1][t2];
		output<<endl;
	}
}

//=====================print=========================

void cell :: print_all(ofstream & output)
{
	// bond length
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_length[t1][t2];
	// coordination number
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_coord_num_surf[t1][t2];
	// composition
	for (int t1=0; t1<num_element; t1++)
		output<<setw(9)<<setprecision(5)<<fixed<<composition_surf[t1];
	// norm of bond vector
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_vec_norm[t1][t2];
	// bond vector in z direction
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_vec_z[t1][t2];
	// global instability index
	output<<setw(9)<<setprecision(5)<<fixed<<gii;
	// global bv
	output<<setw(9)<<setprecision(5)<<fixed<<bv;
	// global bvv
	output<<setw(9)<<setprecision(5)<<fixed<<bvv;
	// total energy
	output<<setw(13)<<setprecision(5)<<fixed<<tot_energy<<endl;
}

void cell :: print_all_sdd(ofstream & output)
{
	// bond length
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_length_sdd[t1][t2];
	// coordination number
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_coord_num_surf_sdd[t1][t2];
	// norm of bond vector
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_vec_norm_sdd[t1][t2];
	// bond vector in z direction
	for (int t1=0; t1<num_element; t1++)
		for (int t2=0; t2<num_element; t2++)
			output<<setw(9)<<setprecision(5)<<fixed<<mean_bond_vec_z_sdd[t1][t2];
	output<<endl;
}
