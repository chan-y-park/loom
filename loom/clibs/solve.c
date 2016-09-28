#include <stdio.h>
#include <complex.h>
#include "solve.h"

newton_params f_df_dx_0(
    diff_params phi,
    double complex z_0,
    double complex x_0
) {
    double complex f_0 = 0;
    double complex df_dx_0 = 0;

    int k;
    double complex c;
    double e;

    newton_params newton;

    for(int i = 0; i < phi.n; i++) {
        k = phi.k[i];
        c = phi.c[i];
        e = phi.e[i];
        f_0 += cpow(x_0, k) * c * cpow(z_0, e);
        if (k > 0) {
            df_dx_0 += k * cpow(x_0, (k - 1)) * c * cpow(z_0, e); 
        }
    }
    newton.f = f_0;
    newton.df_dx = df_dx_0;
    return newton;
}

double complex get_x(
    diff_params phi,
    double complex z_0,
    double complex x_0,
    double accuracy,
    int max_steps
) {
    int step = 0;
    double complex x_i = x_0;
    newton_params newton;
    double complex Delta;

//    printf("x_0 = (%.8f)+(%.8f)I\n", creal(x_0), cimag(x_0));
    while (step < max_steps) {
        newton = f_df_dx_0(phi, z_0, x_i);
        Delta = newton.f / newton.df_dx;
        x_i -= Delta;
/*
        printf(
            "%d:\tf = (%.8f)+(%.8f)I\n"
            "\tdf_dx = (%.8f)+(%.8f)I\n"
            "\tDelta = (%.8f)+(%.8f)I\n",
            step,
            creal(newton.f), cimag(newton.f),
            creal(newton.df_dx), cimag(newton.df_dx),
            creal(Delta), cimag(Delta)
        );
*/
        if (cabs(Delta) < accuracy) {
            break;
        } else {
            step++;
        }
    }
//    printf("x_i = (%.8f)+(%.8f)I\n", creal(x_i), cimag(x_i));
    return x_i;
}
