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
	void clean();			// reset value to zero
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec & operator=(const vec&);
	vec & operator=(double*);	// can replace import
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
	double *composition_surf;
	double **mean_bond_vec_norm;
	double **mean_bond_vec_z;
	double gii;					// global instability index, note that this is not transferrable when changeing the system
	double bv;
	double bvv;
	// sdd of identifer
	double **mean_bond_length_sdd;
	double **mean_coord_num_surf_sdd;
	double **mean_bond_vec_norm_sdd;
	double **mean_bond_vec_z_sdd;
public:
	// end of work
	void clean();
	// 
	void read_file(std::ifstream&);
	void build_e_l(std::string*, int);
	void find_neighbor(double);		//distance between each two elements, could be matrix if necessary
	void find_neighbor(double**);	//distance between each two elements, could be matrix if necessary
	// calculate identifier
	void get_mean_bond_length();
	void get_mean_coord_num_surf(double h_surf);
	void get_composition_surf(double h_surf);
	void get_mean_bond_vec_norm(double h_surf);	// bond vec z included
	void get_gii(double h_surf, std::string el1, std::string el2, double r0, double cc, double n_ox1, double n_ox2);
	void get_bv(double h_surf, std::string el1, std::string el2, double r0, double cc);
	void get_bvv(double h_surf, std::string el1, std::string el2, double r0, double cc);
	// image
	void generate_image(int nx, int ny, std::ofstream&, double *charge, double h_surf);
	// print function
	void print_all(std::ofstream&);
	void print_all_sdd(std::ofstream&);
};

#endif
