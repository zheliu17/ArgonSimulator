#ifndef HEADER_H
#define HEADER_H

using std::string;
using std::vector;

// base.cpp
class Atom
{
public:
    Atom(string atomSymbol,
         vector<double> position,
         vector<double> velocity,
         double mass = 1,
         vector<double> force = {0, 0, 0});
    string symbol;
    vector<double> position;
    vector<double> velocity;
    double mass;
    vector<double> force;
};

class Molecule
{
public:
    Molecule(vector<Atom> atomList = {});
    vector<Atom> atoms;
    void addAtom(Atom atom);
    void printMolecule(std::ostream &output = std::cout) const;
    vector<int> neighbor_list;
    vector<vector<double>> image_position;

    // is_pbc = {}, no pbc;
    // 1, this neighbor mol is not image, regular mol;
    // 2, this neighbor atom is at image_position;
    vector<int> is_pbc;
};

void dumpXYZ(const std::string &XYZfileName,
             const vector<Molecule> &molList,
             bool appendMode = false,
             const std::string &comment = "");

typedef Molecule (*create_mol_func)();
Molecule create_Ar_mol();
void init_mol_list(vector<Molecule> &molList,
                   create_mol_func func,
                   int num_mol = 27,
                   double box_size = 10);

double get_kin_energy(const vector<Molecule> &molList);
double get_temperature(const double &E_kin, const int &natom);
void init_velocity(vector<Molecule> &molList, const double &temperature = 300);

class L_J_potential
{
private:
    double epsilon;
    double sigma;
    double sigma_6;                   // sigma^6
    double sigma_12;                  // sigma^12
    double epsilon_sigma_6_times_24;  // 24 \times \epsilon \times sigma^6
    double epsilon_sigma_12_times_48; // 48 \times \epsilon \times sigma^12
    double epsilon_sigma_6_times_4;   // 4 \times \epsilon \times sigma^6
    double epsilon_sigma_12_times_4;  // 4 \times \epsilon \times sigma^12

public:
    L_J_potential(double i = 1, double j = 1);
    void displayValues();
    double energy(double inv_r);
    double force(double inv_r);
};

// test.cpp
void test_neighbor_list();
void test_neighbor_force();
void test_pbc();
void test_neighbor_search_freq(string filename, const int frequency, const bool is_pbc = false, bool larger = false, int nstep = 2000);
void NVE_half_list(string filename, bool larger = false, int nstep = 2000);
void NVT_grid_search(string filename, const int frequency, const bool is_pbc);
void benchmark_neighbor_freq_pbcon();
void benchmark_neighbor_freq_nopbc();
void test_omp();
void NVT_results();

// inline_function.h
void leap_frog(vector<double> &q,
               vector<double> &v,
               const vector<double> &acc,
               const double dt);

void leap_frog_scale(vector<double> &q,
                     vector<double> &v,
                     const vector<double> &acc,
                     const double dt, const double &lambda);

// misc.cpp
int nearest_cubic(int number);
void sortVectorsByFirst(
    vector<int> &firstVector,
    vector<vector<double>> &secondVector,
    vector<int> &thirdVector);

// neighbor_list.cpp
void clean_neighbor_list(vector<Molecule> &mollist);
void full_neighbor_list(vector<Molecule> &mollist);
void half_neighbor_list(vector<Molecule> &mollist);

class GridSearch
{
public:
    vector<Molecule> &molList;
    double &extend_Rc; // cut-off radius + buffer size
    bool is_pbc_on;    // periodic boundary conditions
    double box_size;
    double Rg; // size of grid, 1/2 extend_Rc
    int numBins;
    vector<vector<vector<vector<int>>>> binIndices;
    double r2;
    double inv_Rg;

    GridSearch(vector<Molecule> &molList, double &R, bool is_pbc_on, double box_size);
    void kernel();

private:
    int getBinIndex(double &position, double inv_Rn, int numBins);
    void processNeighbor(int &neighbor_i, int &neighbor_j, int &neighbor_k, const vector<int> &current_box);
    int getPBCIndex(const int &binIndex, double &off_set);
    void processNeighborsForI(int i, int j, int k, const vector<int> &current_box);
    void processNeighborsForJ(int i, int j, int k, const vector<int> &current_box);
    void processNeighborsForK(int i, int j, int k, const vector<int> &current_box);
    void search_neighbor_box(const vector<int> &current_box, const vector<int> &neighbor_box,
                             const double &off_set_x = 0, const double &off_set_y = 0, const double &off_set_z = 0);
    void clean_binIndices();
};

// algorithm.cpp
class calculate_lj
/* Class to calculate LJ potential force and energy for each mol in a molList
Current only work for single type of mol (water, Ar, etc).
Not work for any mixture*/
{
public:
    vector<Molecule> &molList;
    vector<vector<L_J_potential>> force_matrix;
    /* need to be consistent with mol index
    Example, if H2O is defined as [O, H, H]
    force_matrix[0][0] is the LJ potential between OO */
    bool is_pbc_on;
    size_t natom;
    calculate_lj(vector<Molecule> &molList,
                 const vector<vector<L_J_potential>> &force_matrix,
                 const bool &is_pbc_on = false);

    void kernel(const bool &update_energy = false);
    double get_E_lj();

private:
    double energy;
    void atom_to_atom_LJ(const double &dx, const double &dy, const double &dz,
                         double &fx, double &fy, double &fz,
                         double &temp_fx, double &temp_fy, double &temp_fz,
                         const int &atom_i, const int &atom_j, const bool &update_energy = false);
};

void one_step_md(calculate_lj &calculator, vector<Molecule> &molList, const double &dt, const bool &update_energy);
void one_step_md_scale(calculate_lj &calculator,
                       vector<Molecule> &molList, const double &dt,
                       const bool &update_energy, const double &lambda);
double Berendsen(const double &current_T,
                 const double &target_T,
                 const double &dt_over_tau);

#include "inline_function.h"

#endif
