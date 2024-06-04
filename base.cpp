#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <algorithm>
#include "header.h"

using std::string;
using std::vector;

Atom::Atom(string atomSymbol,
           vector<double> position,
           vector<double> velocity,
           double mass,
           vector<double> force)
    : symbol(atomSymbol), position(position), velocity(velocity), mass(mass), force(force) {}

Molecule::Molecule(vector<Atom> atomList) : atoms(atomList)
{
    is_pbc = {};
}

void Molecule::addAtom(Atom atom)
{
    atoms.push_back(atom);
}

void Molecule::printMolecule(std::ostream &output) const
{
    output << "Molecule Information:" << std::endl;
    for (const Atom &atom : atoms)
    {
        output << "atom Symbol: " << atom.symbol << std::endl;
        output << "Position: [" << atom.position[0] << ", "
               << atom.position[1] << ", " << atom.position[2] << "]" << std::endl;
        output << "Velocity: [" << atom.velocity[0] << ", "
               << atom.velocity[1] << ", " << atom.velocity[2] << "]" << std::endl;
        output << "Force: [" << atom.force[0] << ", "
               << atom.force[1] << ", " << atom.force[2] << "]" << std::endl;
    }
    output << "Neighbor List:" << std::endl;
    for (size_t i = 0; i < neighbor_list.size(); i++)
    {
        output << neighbor_list[i] << " ";
    }
    output << '\n';
    if (!is_pbc.empty())
    {
        for (size_t i = 0; i < neighbor_list.size(); i++)
        {
            if (is_pbc[i] == 1)
            {
                output << "R "; // ith neighbor Mol is a regular Mol
            }
            else
            {
                // ith neighbor Mol is in the image position
                output << "("
                       << image_position[i][0] << ","
                       << image_position[i][1] << ","
                       << image_position[i][2] << ") ";
            }
        }
    }
    output << '\n'
           << "----------------------" << std::endl;
}

void dumpXYZ(const std::string &XYZfileName,
             const vector<Molecule> &molList,
             bool appendMode,
             const std::string &comment)
/* Function to dump molList to XYZ file*/
{
    std::ofstream outputFile;
    if (appendMode)
    {
        outputFile.open(XYZfileName, std::ios::app);
    }
    else
    {
        outputFile.open(XYZfileName);
    }

    outputFile << molList.size() << std::endl;
    outputFile << std::fixed << std::setprecision(6);
    outputFile << comment << std::endl;

    for (const Molecule &mol : molList)
    {
        for (const Atom &atom : mol.atoms)
        {
            outputFile << atom.symbol << "  "
                       << std::setw(10) << atom.position[0] << "  "
                       << std::setw(10) << atom.position[1] << "  "
                       << std::setw(10) << atom.position[2] << std::endl;
        }
    }
    outputFile.close();
}

Molecule create_Ar_mol()
{
    string symbol = "Ar";
    Atom Ar(symbol, {0, 0, 0}, {0, 0, 0}, 39.792, {0, 0, 0});
    Ar.position = {0, 0, 0};
    Molecule Ar_mol(vector<Atom>{Ar});
    return Ar_mol;
}

void init_mol_list(vector<Molecule> &molList,
                   create_mol_func func,
                   int num_mol,
                   double box_size)
/* Function to initialize molList, create evenly distributed molecules
Example, to pop 256 Ar, we created 343 (7*7*7) space point
Then randomly choose 256 of them to fill a atom */
{
    std::mt19937 rng(std::random_device{}());

    int mol_per_axis = nearest_cubic(num_mol);
    double mol_init_dist = box_size / (mol_per_axis);
    double half_mol_init_dist = mol_init_dist / 2;
    vector<bool> is_mol(mol_per_axis * mol_per_axis * mol_per_axis, false);

    for (int i = 0; i < num_mol; ++i)
    {
        is_mol[i] = true;
    }

    std::shuffle(is_mol.begin(), is_mol.end(), rng);

    molList.reserve(num_mol);

    for (int i = 0; i < mol_per_axis; i += 1)
    {
        for (int j = 0; j < mol_per_axis; j += 1)
        {
            for (int k = 0; k < mol_per_axis; k += 1)
            {
                int full_idx = i * mol_per_axis * mol_per_axis + j * mol_per_axis + k;
                if (is_mol[full_idx])
                {
                    Molecule Mol = func();
                    for (size_t l = 0; l < Mol.atoms.size(); l += 1)
                    {
                        Mol.atoms[l].position[0] += half_mol_init_dist + i * mol_init_dist;
                        Mol.atoms[l].position[1] += half_mol_init_dist + j * mol_init_dist;
                        Mol.atoms[l].position[2] += half_mol_init_dist + k * mol_init_dist;
                    }
                    molList.push_back(Mol);
                }
            }
        }
    }
}

double get_kin_energy(const vector<Molecule> &molList)
// Units, kcal per mol
{
    double E_kin = 0;
    for (Molecule mol : molList)
    {
        for (Atom atom : mol.atoms)
        {
            vector<double> v = atom.velocity;
            E_kin += (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * atom.mass * 0.5;
        }
    }
    E_kin = E_kin * 1e4 / 4.184; // from g * (A per fs)^2 to kcal
    return E_kin;
}

double get_temperature(const double &E_kin, const int &natom)
// Units, K 
// By default, consider 3n - 6 degree of freedom
{
    double kb = 1.987204259e-3; // Boltzmann constant kcal per mol
    return 2 * E_kin / (3 * natom - 6) / kb;
}

void init_velocity(vector<Molecule> &molList, const double &temperature)
// initialize velocity according to Maxwell-Boltzmann distribution
{
    std::random_device entropy_source;
    std::mt19937_64 gen(entropy_source());
    double R = 8.3145; // Gas constant
    // Consider 3n - 6 degree of freedom
    double Kb_T = R * temperature / 1e7 * (1 + 2 / molList.size());
    for (Molecule &mol : molList)
    {
        for (Atom &atom : mol.atoms)
        {
            double stddev = sqrt(Kb_T / atom.mass);
            std::normal_distribution<double> dist(0, stddev);
            atom.velocity = {dist(gen), dist(gen), dist(gen)};
        }
    }
}

L_J_potential::L_J_potential(double i, double j) : epsilon(i), sigma(j)
{
    sigma_6 = pow(sigma, 6);
    sigma_12 = sigma_6 * sigma_6;
    epsilon_sigma_6_times_24 = epsilon * sigma_6 * 24;
    epsilon_sigma_12_times_48 = epsilon * sigma_12 * 48;
    epsilon_sigma_6_times_4 = epsilon * sigma_6 * 4;
    epsilon_sigma_12_times_4 = epsilon * sigma_12 * 4;
}

void L_J_potential::displayValues()
{
    std::cout << "epsilon: " << epsilon << ", sigma: " << sigma << std::endl;
}
