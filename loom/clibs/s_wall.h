typedef struct {int s_wall_size; int rv;} message;
typedef struct {double complex _1; double complex _2;} ode_xs; 
typedef struct {
    double size_of_small_step;
    double size_of_large_step;
    double size_of_bp_neighborhood;
    double size_of_puncture_cutoff;
    double mass_limit;
    double accuracy;
} numerical_parameters;
typedef struct {double start; double end;} twist_line;

#define ERROR_SAME_XS -1;
