#include <stdio.h>
#include <complex.h>
#include "s_wall.h"
#include "solve.h"

int grow(
    message* msg,
    int N,
    int* diff_n_k, double complex* diff_n_c, double* diff_n_e, int n_diff_n,
    int* diff_d_k, double complex* diff_d_c, double* diff_d_e, int n_diff_d,
    double complex* c_dz_dt,
    double complex* z,
    ode_xs* x,
    double* M,
    double complex* bpz, int n_bpz,
    double complex* ppz, int n_ppz,
    numerical_parameters np,
    twist_line* tl, int n_tl
) {
    int s_wall_size = msg->s_wall_size;
    int max_steps = msg->rv;

    diff_params phi_n;
    phi_n.k = diff_n_k;
    phi_n.c = diff_n_c;
    phi_n.e = diff_n_e;
    phi_n.n = n_diff_n;

    diff_params phi_d;
    phi_d.k = diff_d_k;
    phi_d.c = diff_d_c;
    phi_d.e = diff_d_e;
    phi_d.n = n_diff_d;

    double complex z_i, x_i_1, x_i_2, x_tmp;
    double M_i;

    int i = 0;
    double min_d;
    double d;
    double dt;
    double complex Dx_i;
    double avg_z_r;

    while (i < (s_wall_size - 1)) {
        z_i = z[i];
        x_i_1 = x[i]._1;
        x_i_2 = x[i]._2;
        M_i = M[i];

        i++;
        if (i > MIN_NUM_OF_DATA_PTS) {
            min_d = np.size_of_puncture_cutoff;
            for (int j = 0; j < n_ppz; j++) {
                d = cabs(z_i - ppz[j]);
                if (d < min_d) min_d = d;
            }
            if (min_d < np.size_of_puncture_cutoff) {
                msg->s_wall_size = i;
                msg->rv = NEAR_PUNCTURE;
                return 0;
            }

            if (M_i > np.mass_limit) {
                msg->s_wall_size = i;
                msg->rv = MASS_LIMIT;
                return 0;
            }
        }

        min_d = np.size_of_bp_neighborhood;
        for (int j = 0; j < n_bpz; j++) {
            d = cabs(z_i - bpz[j]);
            if (d < min_d) min_d = d;
        }
        if (min_d < np.size_of_bp_neighborhood) {
            dt = np.size_of_small_step;
        } else {
            dt = np.size_of_large_step;
        }

        Dx_i = x_i_1 - x_i_2;
        //z[i] = z_i + dt * c_dz_dt[0] / Dx_i;
        z[i] = z_i + dt * cexp((carg(c_dz_dt[0]) - carg(Dx_i))*I);
        M[i] = M_i + cabs(Dx_i) * dt;

        if (n_tl > 0 && cimag(z[i]) * cimag(z_i) < 0) {
            avg_z_r = (creal(z[i]) + creal(z_i)) * 0.5;
            for (int j = 0; j < n_tl; j++) {
                if (tl[j].start < avg_z_r && avg_z_r < tl[j].end) {
                    x_tmp = x_i_1;
                    x_i_1 = -1.0 * x_i_2;
                    x_i_2 = -1.0 * x_tmp;
                    break;
                }
            }
        }

        x[i]._1 = get_x(N, phi_n, phi_d, z[i], x_i_1, np.accuracy, max_steps);
        x[i]._2 = get_x(N, phi_n, phi_d, z[i], x_i_2, np.accuracy, max_steps);

        if (cabs(x[i]._1 - x[i]._2) < np.accuracy) {
            msg->s_wall_size = i;
            msg->rv = ERROR_SAME_XS;
            return 0;
        }

    }

    return 0;
}
