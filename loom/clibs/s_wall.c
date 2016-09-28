#include <stdio.h>
#include <complex.h>
#include "s_wall.h"
#include "solve.h"

int grow(
    message* msg,
    int* diff_k, double complex* diff_c, double* diff_e, int n_diff,
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

    diff_params phi;
    phi.k = diff_k;
    phi.c = diff_c;
    phi.e = diff_e;
    phi.n = n_diff;

    double complex z_i, x_i_1, x_i_2, x_tmp;
    double M_i;
    //double complex x_n_1, x_n_2;

/*
    int i;
    for (i = 0; i < n_diff; i++) {
        k = diff_k[i];
        c = diff_c[i];
        e = diff_e[i];
        printf(
            "(k, c, e) = (%d, (%.8f)+(%.8f)I, %f)\n",
            k, creal(c), cimag(c), e
        );
    }

    printf(
        "c_dz_dt[0] = (%.8f)+(%.8f)I\n",
        creal(c_dz_dt[0]), cimag(c_dz_dt[0])
    );

    printf("z[0] = (%.8f)+(%.8f)I\n", creal(z[0]), cimag(z[0]));
    printf("x[0][0] = (%.8f)+(%.8f)I\n", creal(x[0]._1), cimag(x[0]._2));
    printf("x[0][1] = (%.8f)+(%.8f)I\n", creal(x[0]._2), cimag(x[0]._2));
    printf("M[0] = (%.8f)+(%.8f)I\n", creal(M[0]), cimag(M[0]));

    for (i = 0; i < n_bpz; i++) {
        printf("bpz[%d] = (%.8f)+(%.8f)I\n", i, creal(bpz[i]), cimag(bpz[i]));
    }

    for (i = 0; i < n_ppz; i++) {
        printf("ppz[%d] = (%.8f)+(%.8f)I\n", i, creal(ppz[i]), cimag(ppz[i]));
    }

    printf("size_of_small_step = %.8f\n", np.size_of_small_step);
    printf("size_of_large_step = %.8f\n", np.size_of_large_step);
    printf("size_of_bp_neighborhood = %.8f\n", np.size_of_bp_neighborhood);
    printf("size_of_puncture_cutoff = %.8f\n", np.size_of_puncture_cutoff);
    printf("mass_limit = %.8f\n", np.mass_limit);
    printf("accuracy = %.8f\n", np.accuracy);
    printf("s_wall_size = %d\n", msg->s_wall_size);
    printf("rv = %d\n", msg->rv);
*/

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
/*
        printf(
            "%d: M[i] = %.8f, M_i = %.8f, dt = %.8f\n",
            i, M[i], M_i, dt
        );
*/
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

        x[i]._1 = get_x(phi, z[i], x_i_1, np.accuracy, max_steps);
        x[i]._2 = get_x(phi, z[i], x_i_2, np.accuracy, max_steps);

        if (cabs(x[i]._1 - x[i]._2) < np.accuracy) {
            msg->s_wall_size = i;
            msg->rv = ERROR_SAME_XS;
            return 0;
        }

    }

    return 0;
}
