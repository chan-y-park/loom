typedef struct {int* k; double complex* c; double* e; int n;} diff_params;
typedef struct {double complex f; double complex df_dx;} newton_params;

newton_params f_df_dx_0(diff_params, double complex, double complex);

double get_x(diff_params, double complex, double complex, double, int);
