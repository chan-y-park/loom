# Factors regarding various sizes.
# F_SIZE_OF_SMALL_STEPS = 1e-3
F_SIZE_OF_SMALL_STEPS = 1e-2
# F_SIZE_OF_LARGE_STEPS = 1e-2
F_SIZE_OF_LARGE_STEPS = 1e-1
F_SIZE_OF_BP_NEIGHBORHOOD = .5
F_SIZE_OF_PUNCTURE_CUTOFF = 1e-2

# Large numbers for infinity.
N_INF = -1e308
P_INF = +1e308

# S-wall grow method libraries.
LIB_C = 0
LIB_SCIPY_ODE = 1
LIB_NUMBA = 2
LIB_PYTHON = 3

# NOTE: The following should have the same values as those in clibs/s_wall.h
ERROR_SAME_XS = -1
NEAR_PUNCTURE = 1
MASS_LIMIT = 2
IN_BP_NBHD = 3
OUT_BP_NBHD = 4

# Additional error code for SciPy ODE
ERROR_SCIPY_ODE = -2

# Maximum # of steps when using Newton's method.
NEWTON_MAX_STEPS = 100

# Default min/max # of data points per S-wall.
S_WALL_MIN_N_DATA_PTS = 100
S_WALL_MAX_N_DATA_PTS = 200

# Number of iterations when improving a soliton tree
SOLITON_TREE_MAX_N_ITERS = 10

# Minimum number of dz steps when improving a soliton tree
SOLITON_TREE_MIN_N_DZ_STEPS = 10
