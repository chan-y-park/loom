typedef struct {int* k; double complex* c; double* e; int n;} diff_params;
typedef struct {double complex f; double complex df_dx;} newton_params;

newton_params f_df_dx_0(
    int, diff_params, diff_params, double complex, double complex
);

double complex get_x(
    int, diff_params, diff_params, double complex, double complex, double, int
);
