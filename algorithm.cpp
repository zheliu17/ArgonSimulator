#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "header.h"

using std::vector;

calculate_lj::calculate_lj(vector<Molecule> &molList,
                           const vector<vector<L_J_potential>> &force_matrix,
                           const bool &is_pbc_on)
    : molList(molList),
      force_matrix(force_matrix),
      is_pbc_on(is_pbc_on),
      natom(molList[0].atoms.size()),
      energy(0) {}

void calculate_lj::kernel(const bool &update_energy)
{
    if (update_energy)
    {
        energy = 0;
    }
    // clear force
#pragma omp parallel for
    for (Molecule &mol : molList)
    {
        for (Atom &atom : mol.atoms)
        {
            atom.force = {0, 0, 0};
        }
    }

#pragma omp parallel for
    for (size_t mol_i = 0; mol_i < molList.size(); mol_i++)
    {
        const vector<int> &neighbor_list_ref = molList[mol_i].neighbor_list;
        const size_t &neighbor_list_size = neighbor_list_ref.size();

        for (size_t atom_j = 0; atom_j < natom; atom_j++)
        {
            Atom &atom_j_ref = molList[mol_i].atoms[atom_j];
            double x0 = atom_j_ref.position[0];
            double y0 = atom_j_ref.position[1];
            double z0 = atom_j_ref.position[2];
            double temp_fx = 0;
            double temp_fy = 0;
            double temp_fz = 0;
            for (size_t k = 0; k < neighbor_list_size; k++)
            {
                for (size_t atom_l = 0; atom_l < natom; atom_l++)
                {
                    Atom &atom_l_ref = molList[neighbor_list_ref[k]].atoms[atom_l];
                    double dx = atom_l_ref.position[0] - x0;
                    double dy = atom_l_ref.position[1] - y0;
                    double dz = atom_l_ref.position[2] - z0;
                    if (is_pbc_on && molList[mol_i].is_pbc[k] == 2)
                    {
                        dx += molList[mol_i].image_position[k][0];
                        dy += molList[mol_i].image_position[k][1];
                        dz += molList[mol_i].image_position[k][2];
                    }
                    double fx = 0;
                    double fy = 0;
                    double fz = 0;
                    atom_to_atom_LJ(dx, dy, dz,
                                    fx, fy, fz,
                                    temp_fx, temp_fy, temp_fz,
                                    atom_j, atom_l, update_energy);
#pragma omp atomic
                    atom_l_ref.force[0] -= fx;
#pragma omp atomic
                    atom_l_ref.force[1] -= fy;
#pragma omp atomic
                    atom_l_ref.force[2] -= fz;
                }
            }
#pragma omp atomic
            atom_j_ref.force[0] += temp_fx;
#pragma omp atomic
            atom_j_ref.force[1] += temp_fy;
#pragma omp atomic
            atom_j_ref.force[2] += temp_fz;
        }
    }
}

double calculate_lj::get_E_lj()
{
    return energy;
}

void calculate_lj::atom_to_atom_LJ(const double &dx, const double &dy, const double &dz,
                                   double &fx, double &fy, double &fz,
                                   double &temp_fx, double &temp_fy, double &temp_fz,
                                   const int &atom_i, const int &atom_j, const bool &update_energy)
{
    double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    double inv_r = 1 / r;
    double temp_f = force_matrix[atom_i][atom_j].force(inv_r) * inv_r;
    if (update_energy)
    {
#pragma omp atomic
        energy += force_matrix[atom_i][atom_j].energy(inv_r);
    }
    fx = temp_f * dx;
    fy = temp_f * dy;
    fz = temp_f * dz;
    temp_fx += fx;
    temp_fy += fy;
    temp_fz += fz;
}

void one_step_md(calculate_lj &calculator, vector<Molecule> &molList, const double &dt, const bool &update_energy)
{
    calculator.kernel(update_energy);
    double kcal_per_mol_per_A_to_A_per_fs = 4184 / 1e7;
#pragma omp parallel for
    for (size_t i = 0; i < molList.size(); i++)
    {
        for (Atom &atom : molList[i].atoms)
        {
            vector<double> acc(3);
            double factor = kcal_per_mol_per_A_to_A_per_fs / atom.mass;
            acc[0] = atom.force[0] * factor;
            acc[1] = atom.force[1] * factor;
            acc[2] = atom.force[2] * factor;
            leap_frog(atom.position, atom.velocity, acc, dt);
        }
    }
}

void one_step_md_scale(calculate_lj &calculator,
                       vector<Molecule> &molList, const double &dt,
                       const bool &update_energy, const double &lambda)
{
    calculator.kernel(update_energy);
    double kcal_per_mol_per_A_to_A_per_fs = 4184 / 1e7;

#pragma omp parallel for
    for (size_t i = 0; i < molList.size(); i++)
    {
        for (Atom &atom : molList[i].atoms)
        {
            vector<double> acc(3);
            double factor = kcal_per_mol_per_A_to_A_per_fs / atom.mass;
            acc[0] = atom.force[0] * factor;
            acc[1] = atom.force[1] * factor;
            acc[2] = atom.force[2] * factor;
            leap_frog_scale(atom.position, atom.velocity, acc, dt, lambda);
        }
    }
}

double Berendsen(const double &current_T,
                 const double &target_T,
                 const double &dt_over_tau)
// Berendsen thermostat
{
    return sqrt(1 + dt_over_tau * (target_T / current_T - 1));
}