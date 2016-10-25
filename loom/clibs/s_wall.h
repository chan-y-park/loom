typedef struct {
    int s_wall_size;
    int step;
    int stop_condition;
    int rv;
} message;
typedef struct {double complex _1; double complex _2;} ode_xs; 
typedef struct {
    double size_of_small_step;
    double size_of_large_step;
    double size_of_bp_neighborhood;
    double size_of_pp_neighborhood;
    double size_of_puncture_cutoff;
    double mass_limit;
    double accuracy;
} numerical_parameters;
typedef struct {double start; double end;} twist_line;

#define ERROR_SAME_XS -1
#define NEAR_PUNCTURE 1
#define MASS_LIMIT 2
#define IN_P_NBHD 3
#define OUT_P_NBHD 4

#define MIN_NUM_OF_DATA_PTS 3

#define SAME_XS_MAX_STEPS 5

#define N_INF -1e308
#define P_INF +1e308

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
);


double get_min_d(double complex z, double complex* pz, int n_pz);
