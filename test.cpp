#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <algorithm>
#include "header.h"
#include <omp.h>

using std::string;
using std::vector;
using std::chrono::duration;
using std::chrono::high_resolution_clock;

void test_neighbor_list()
// Please refers to NVT_grid_search() for the choice of the constants
{
    vector<Molecule> molList;
    double box_size = 22.9;
    double Rc = 8;
    double delta_s = 2;
    double extend_Rc = Rc + delta_s;
    std::ofstream test("debug.log");
    std::ofstream debug_output("debug.log", std::ios::app);

    init_mol_list(molList, create_Ar_mol, 256, box_size);
    string filename = "out.xyz";
    dumpXYZ(filename, molList, false);

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    full_neighbor_list(molList);

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    half_neighbor_list(molList);
    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    GridSearch GridSearcher(molList, extend_Rc, true, box_size);
    GridSearcher.kernel();

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }
    std::cout << "debug file -> \"debug.log\"" << std::endl;
}

void test_neighbor_force()
// Please refers to NVT_grid_search() for the choice of the constants
{
    double epsilon = 1.987204259e-3 * 119.8;
    L_J_potential RNP(epsilon, 3.405);

    vector<Molecule> molList;
    double box_size = 22.9;
    double Rc = 8;
    double delta_s = 2;
    double extend_Rc = Rc + delta_s;

    std::ofstream test("debug.log");
    std::ofstream test2("debug2.log");
    std::ofstream debug_output("debug.log", std::ios::app);
    std::ofstream debug_output2("debug2.log", std::ios::app);

    init_mol_list(molList, create_Ar_mol, 256, box_size);
    string filename = "out.xyz";
    dumpXYZ(filename, molList, false);

    half_neighbor_list(molList);
    vector<L_J_potential> fm1 = {RNP};
    vector<vector<L_J_potential>> full_fm = {fm1};
    calculate_lj calculator(molList, full_fm, false);
    calculator.kernel();

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    GridSearch GridSearcher(molList, extend_Rc, false, box_size);
    GridSearcher.kernel();
    calculator.kernel();

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output2 << "index " << i << std::endl;
        molList[i].printMolecule(debug_output2);
    }
    std::cout << "debug file -> \"debug.log\" and \"debug2.log\"" << std::endl;
};

void test_pbc()
// Please refers to NVT_grid_search() for the choice of the constants
{
    double epsilon = 1.987204259e-3 * 119.8;
    L_J_potential RNP(epsilon, 3.405);
    vector<vector<L_J_potential>> full_fm = {{RNP}};

    vector<Molecule> molList;
    double box_size = 22.9;
    Molecule mol = create_Ar_mol();
    mol.atoms[0].position = {1, 1, 1};
    molList.push_back(mol);
    mol.atoms[0].position = {11, 11, 11};
    molList.push_back(mol);
    mol.atoms[0].position = {21.9, 1, 1};
    molList.push_back(mol);
    double Rc = 8;
    double delta_s = 2;
    double extend_Rc = Rc + delta_s;

    std::ofstream test("debug.log");
    std::ofstream debug_output("debug.log", std::ios::app);

    string filename = "out.xyz";
    dumpXYZ(filename, molList, false);

    GridSearch GridSearcher(molList, extend_Rc, false, box_size);
    GridSearcher.kernel();

    calculate_lj calculator(molList, full_fm, false);
    calculator.kernel();

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    GridSearcher.is_pbc_on = true;
    GridSearcher.kernel();
    calculator.is_pbc_on = true;
    calculator.kernel();

    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    molList[2].atoms[0].position = {-1, 1, 1};
    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }

    GridSearcher.kernel();
    calculator.kernel();
    for (size_t i = 0; i < molList.size(); i += 1)
    {
        debug_output << "index " << i << std::endl;
        molList[i].printMolecule(debug_output);
    }
}

void NVE_half_list(string filename, bool larger, int nstep)
// This option is to use half_neighbor_list as the neighbor list and never update the list. 
// Which essentially means calculating force between all atoms without any cut-off 
// Please refers to NVT_grid_search() for the choice of the constants
{
    double epsilon = 1.987204259e-3 * 119.8;
    L_J_potential RNP(epsilon, 3.405);
    vector<vector<L_J_potential>> full_fm = {{RNP}};

    vector<Molecule> molList;
    double box_size = 22.9;
    if (larger)
    {
        box_size = 22.9 * 2;
    }
    double dt = 0.01; // fs

    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    duration<double, std::milli> duration_sec;

    std::ofstream file1("out.xyz", std::ios::trunc);
    std::ofstream file2(filename, std::ios::trunc);
    std::ofstream debug_output(filename, std::ios::app);

    int natom = 256;
    if (larger)
    {
        natom = 256 * 8;
    }
    init_mol_list(molList, create_Ar_mol, natom, box_size);
    init_velocity(molList, 85);

    bool is_pbc = false;
    calculate_lj calculator(molList, full_fm, is_pbc);

    double E_kin;

    half_neighbor_list(molList);
    debug_output << " Time (fs)    Kinetic E    Temperature     Potential       Total E" << std::endl;
    debug_output << std::fixed << std::setprecision(6);

    start = high_resolution_clock::now();
    for (int i = 0; i <= nstep; i++)
    {
        if (i % 1000 == 0)
        {
            one_step_md(calculator, molList, dt, true);
            E_kin = get_kin_energy(molList);
            dumpXYZ("out.xyz", molList, true);
            debug_output << std::setw(10) << i * dt << "    "
                         << std::setw(10) << E_kin << "    "
                         << std::setw(10) << get_temperature(E_kin, molList.size()) << "    "
                         << std::setw(10) << calculator.get_E_lj() << "    "
                         << std::setw(10) << E_kin + calculator.get_E_lj() << std::endl;
        }
        else
        {
            one_step_md(calculator, molList, dt, false);
        }
    }
    end = high_resolution_clock::now();

    duration_sec = std::chrono::duration_cast<duration<double, std::milli>>(end - start);
    std::cout << "Elapsed time " << duration_sec.count() << "\n";
}

void test_neighbor_search_freq(string filename, const int frequency, const bool is_pbc, bool larger, int nstep)
// Please refers to NVT_grid_search() for the choice of the constants
{
    double epsilon = 1.987204259e-3 * 119.8;
    L_J_potential RNP(epsilon, 3.405);
    vector<vector<L_J_potential>> full_fm = {{RNP}};

    vector<Molecule> molList;
    double box_size = 22.9;
    if (larger)
    {
        box_size = 22.9 * 2;
    }
    double Rc = 14;
    double delta_s = 2;
    double extend_Rc = Rc + delta_s;
    double dt = 0.01; // fs

    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    duration<double, std::milli> duration_sec;

    std::ofstream file1("out.xyz", std::ios::trunc);
    std::ofstream file2(filename, std::ios::trunc);
    std::ofstream debug_output(filename, std::ios::app);

    int natom = 256;
    if (larger)
    {
        natom = 256 * 8;
    }
    init_mol_list(molList, create_Ar_mol, natom, box_size);
    init_velocity(molList, 85);

    GridSearch GridSearcher(molList, extend_Rc, is_pbc, box_size);
    calculate_lj calculator(molList, full_fm, is_pbc);

    double E_kin;
    debug_output << " Time (fs)    Kinetic E    Temperature     Potential       Total E" << std::endl;
    debug_output << std::fixed << std::setprecision(6);

    start = high_resolution_clock::now();
    for (int i = 0; i <= nstep; i++)
    {
        if (i % frequency == 0)
        {
            GridSearcher.kernel();
        }
        if (i % 1000 == 0)
        {
            one_step_md(calculator, molList, dt, true);
            E_kin = get_kin_energy(molList);
            dumpXYZ("out.xyz", molList, true);
            debug_output << std::setw(10) << i * dt << "    "
                         << std::setw(10) << E_kin << "    "
                         << std::setw(10) << get_temperature(E_kin, molList.size()) << "    "
                         << std::setw(10) << calculator.get_E_lj() << "    "
                         << std::setw(10) << E_kin + calculator.get_E_lj() << std::endl;
        }
        else
        {
            one_step_md(calculator, molList, dt, false);
        }
    }
    end = high_resolution_clock::now();

    duration_sec = std::chrono::duration_cast<duration<double, std::milli>>(end - start);
    std::cout << "Elapsed time " << duration_sec.count() << "\n";
}

void NVT_grid_search(string filename, const int frequency, const bool is_pbc)
// Function to generate the final data
{
    double epsilon = 1.987204259e-3 * 119.8; // kcal per mol
    L_J_potential RNP(epsilon, 3.405);       // J. Comput. Phys. 17.4 (1975): 401-414.
    vector<vector<L_J_potential>> full_fm = {{RNP}};

    vector<Molecule> molList;
    double box_size = 22.9; // Angstrom (10^-10 m)
    // 256 argon at a box size 22.9 correspond to liquid density 1.409 g per cm^3 (Experimental data)
    double Rc = 14;     // cut-off radius
    double delta_s = 2; // buffer size
    double extend_Rc = Rc + delta_s;
    // LJ potential is a short range force, we can truncate it within a cut-off radius
    // add a buffer so that we don't need to update neighbor list every MD step
    double dt = 0.01; // fs

    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    duration<double, std::milli> duration_sec;

    std::ofstream file1("out.xyz", std::ios::trunc); // clean "out.xyz" contains, "out.xyz" will be the trajectory file
    std::ofstream file2(filename, std::ios::trunc);  // dump log to "filename"
    std::ofstream debug_output(filename, std::ios::app);

    init_mol_list(molList, create_Ar_mol, 256, box_size);
    init_velocity(molList, 85); // initialize velocity

    GridSearch GridSearcher(molList, extend_Rc, is_pbc, box_size); // use grid search, and pbc
    calculate_lj calculator(molList, full_fm, is_pbc);

    double E_kin;
    debug_output << " Time (fs)    Kinetic E    Temperature     Potential       Total E" << std::endl;
    debug_output << std::fixed << std::setprecision(6);

    start = high_resolution_clock::now();
    for (int i = 0; i <= 2000000; i++)
    {
        // update neighbor list every "frequency" step
        if (i % frequency == 0)
        {
            GridSearcher.kernel();
        }
        E_kin = get_kin_energy(molList);
        double current_T = get_temperature(E_kin, molList.size());
        double lambda = 1;
        // temperature control
        lambda = Berendsen(current_T, 85, 0.0001);
        // generate log every 1000 step
        if (i % 1000 == 0)
        {
            one_step_md_scale(calculator, molList, dt, true, lambda);

            dumpXYZ("out.xyz", molList, true);
            debug_output << std::setw(10) << i * dt << "    "
                         << std::setw(10) << E_kin << "    "
                         << std::setw(10) << current_T << "    "
                         << std::setw(10) << calculator.get_E_lj() << "    "
                         << std::setw(10) << E_kin + calculator.get_E_lj() << std::endl;
        }
        else
        {
            one_step_md_scale(calculator, molList, dt, false, lambda);
        }
    }
    end = high_resolution_clock::now();

    duration_sec = std::chrono::duration_cast<duration<double, std::milli>>(end - start);
    std::cout << "Elapsed time " << duration_sec.count() << "\n";
}

void benchmark_neighbor_freq_nopbc()
{
    NVE_half_list("debug1.out", true, 2000);
    test_neighbor_search_freq("debug2.out", 1000, false, true, 2000);
    test_neighbor_search_freq("debug3.out", 100, false, true, 2000);
    test_neighbor_search_freq("debug4.out", 10, false, true, 2000);
}

void benchmark_neighbor_freq_pbcon()
{
    test_neighbor_search_freq("debug2.out", 10, true);
    test_neighbor_search_freq("debug3.out", 100, true);
    test_neighbor_search_freq("debug4.out", 1000, true);
}

void test_omp()
{
    bool larger = false;
    omp_set_num_threads(1);
    test_neighbor_search_freq("debug1.out", 100, larger);
    NVE_half_list("debug2.out", larger);
    omp_set_num_threads(4);
    test_neighbor_search_freq("debug3.out", 100, larger);
    NVE_half_list("debug4.out", larger);
    omp_set_num_threads(8);
    test_neighbor_search_freq("debug1.out", 100, larger);
    NVE_half_list("debug2.out", larger);
    omp_set_num_threads(16);
    test_neighbor_search_freq("debug3.out", 100, larger);
    NVE_half_list("debug4.out", larger);

    larger = true;
    omp_set_num_threads(1);
    test_neighbor_search_freq("debug1.out", 100, larger);
    NVE_half_list("debug2.out", larger);
    omp_set_num_threads(4);
    test_neighbor_search_freq("debug3.out", 100, larger);
    NVE_half_list("debug4.out", larger);
    omp_set_num_threads(8);
    test_neighbor_search_freq("debug1.out", 100, larger);
    NVE_half_list("debug2.out", larger);
    omp_set_num_threads(16);
    test_neighbor_search_freq("debug3.out", 100, larger);
    NVE_half_list("debug4.out", larger);
}

void NVT_results()
{
    omp_set_num_threads(4);
    NVT_grid_search("Trajectory.log", 100, true);
}