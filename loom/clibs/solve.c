#include <stdio.h>
#include <complex.h>
#include "solve.h"

newton_params f_df_dx_0(
    int N,
    diff_params phi_n,
    diff_params phi_d,
    double complex z_0,
    double complex x_0
) {
    double complex phi_[N][2];
    double complex f_0 = cpow(x_0, N);
    double complex df_dx_0 = N * cpow(x_0, N - 1);

    int i;
    int k;
    double complex c;
    double e;

    double complex phi_k;

    newton_params newton;

    for(i = 0; i < N; i++) {
        phi_[i][0] = 0;
        phi_[i][1] = 0;
    }

    for(i = 0; i < phi_n.n; i++) {
        k = phi_n.k[i];
        c = phi_n.c[i];
        e = phi_n.e[i];
        phi_[k][0] += c * cpow(z_0, e);
    }

    for(i = 0; i < phi_d.n; i++) {
        k = phi_d.k[i];
        c = phi_d.c[i];
        e = phi_d.e[i];
        phi_[k][1] += c * cpow(z_0, e);
    }

    for(k = 0; k < N; k++) {
        phi_k = phi_[k][0] / phi_[k][1];
        f_0 += phi_k * cpow(x_0, k);
        if (k > 0) {
            df_dx_0 += k * phi_k * cpow(x_0, (k - 1)); 
        }
    }
    newton.f = f_0;
    newton.df_dx = df_dx_0;
    return newton;
}

double complex get_x(
    int N,
    diff_params phi_n,
    diff_params phi_d,
    double complex z_0,
    double complex x_0,
    double accuracy,
    int max_steps
) {
    int step = 0;
    double complex x_i = x_0;
    newton_params newton;
    double complex Delta;

    while (step < max_steps) {
        newton = f_df_dx_0(N, phi_n, phi_d, z_0, x_i);
        Delta = newton.f / newton.df_dx;
        x_i -= Delta;

        if (cabs(Delta) < accuracy) {
            break;
        } else {
            step++;
        }
    }
    return x_i;
}
