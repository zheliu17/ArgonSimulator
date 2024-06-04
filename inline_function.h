#include "header.h"

inline double L_J_potential::energy(double inv_r)
{
    double inv_r_6 = inv_r * inv_r * inv_r * inv_r * inv_r * inv_r;
    double inv_r_12 = inv_r_6 * inv_r_6;
    return inv_r_12 * epsilon_sigma_12_times_4 - inv_r_6 * epsilon_sigma_6_times_4;
}

inline double L_J_potential::force(double inv_r)
// opposite direction of the gradient
{
    double inv_r_6 = inv_r * inv_r * inv_r * inv_r * inv_r * inv_r;
    double inv_r_7 = inv_r_6 * inv_r;
    double inv_r_13 = inv_r_6 * inv_r_7;
    return inv_r_7 * epsilon_sigma_6_times_24 - inv_r_13 * epsilon_sigma_12_times_48;
}

inline void leap_frog(vector<double> &q,
                      vector<double> &v,
                      const vector<double> &acc,
                      const double dt)
{
    v[0] += acc[0] * dt;
    v[1] += acc[1] * dt;
    v[2] += acc[2] * dt;
    q[0] += dt * v[0];
    q[1] += dt * v[1];
    q[2] += dt * v[2];
}

inline void leap_frog_scale(vector<double> &q,
                            vector<double> &v,
                            const vector<double> &acc,
                            const double dt, const double &lambda)
{
    // p-T-p-X
    // first half step update v
    double half_dt = dt * 0.5;
    v[0] += acc[0] * half_dt;
    v[1] += acc[1] * half_dt;
    v[2] += acc[2] * half_dt;
    // then temperature control
    v[0] *= lambda;
    v[1] *= lambda;
    v[2] *= lambda;
    // then second half step
    v[0] += acc[0] * half_dt;
    v[1] += acc[1] * half_dt;
    v[2] += acc[2] * half_dt;
    q[0] += dt * v[0];
    q[1] += dt * v[1];
    q[2] += dt * v[2];
}