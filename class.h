#include <cstring>
#include <fstream>
#ifndef MY_CLASS
#define MY_CLASS

class vec;
class atom;
class cell;

class vec
{
private:
	double x[3];
friend atom;
friend cell;
public:
	void import(double *);
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec & operator=(const vec&);
	vec & operator=(double*);
	double norm();
	//debug
	void print();
};

class atom
{
private:
	std::string symbol;
	int sym;
	vec pos;
friend cell;
public:
	atom & operator=(const atom&);
};

class cell
{
private:
	// atom related
	int num_atom;
	int num_element;
	atom * a_l;			// atom list
	std::string * e_l;		// element list
	int * num_neighbor;
	int max_num_neighbor = 50;
	atom ** neighbor_l;	// neighbor list
	double ** neighbor_dis;
	// cell related
	vec lattice[3];
	double tot_energy;
	// identifier
	double **mean_bond_length;	//between different element, can add height constrain
	double **mean_coord_num_surf;

public:
	// end of work
	void clean();
	// 
	void read_file(std::ifstream&);
	void find_neighbor(double);		//distance between each two elements, could be matrix if necessary
	void find_neighbor(double**);	//distance between each two elements, could be matrix if necessary
	// calculate identifier
	// ---bond length
	void get_mean_bond_length();
	void get_mean_coord_num_surf(double h_surf);
};

#endif
