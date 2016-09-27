#include <complex.h>
#include "solve.h"


newton_params f_df_dx_0(
    diff_params phi,
    double complex x_0,
    double complex z_0
) {
    double f_0 = 0;
    double df_dx_0 = 0;
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
    double complex x_0,
    double complex z_0,
    double accuracy,
    int max_steps
) {
    int step = 0;
    double complex x_i = x_0;
    newton_params newton;
    double complex Delta;

    while (step < max_steps) {
        newton = f_df_dx_0(phi, x_i, z_0);
        Delta = newton.f / newton.df_dx;
        if (cabs(Delta) > accuracy) {
            break;
        }
        x_i -= Delta;
        step++;
    }
    return x_i;
}
