#include <stdio.h>
#include <complex.h>
#include "s_wall.h"

int grow(
    message* msg,
    int* diff_k, double complex* diff_c, double* diff_e, int n_diff,
    double complex* c_dz_dt,
    double complex* z,
    ode_xs* x,
    double complex* M,
    double complex* bpz, int n_bpz,
    double complex* ppz, int n_ppz,
//    double size_of_small_step,
//    double size_of_large_step,
//    double size_of_bp_neighborhood,
//    double size_of_puncture_cutoff,
//    double mass_limit,
//    double accuracy,
    numerical_parameters np,
    twist_line* tl, int n_tl
) {
    int i;
    int k;
    double complex c;
    double e;

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
    printf("x[0][0] = (%.8f)+(%.8f)I\n", creal(x[0].x1), cimag(x[0].x2));
    printf("x[0][1] = (%.8f)+(%.8f)I\n", creal(x[0].x2), cimag(x[0].x2));
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

    msg->s_wall_size -= 10;
    msg->rv = ERROR_SAME_XS;
    return 0;
}
