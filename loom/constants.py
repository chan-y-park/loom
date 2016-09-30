# Factors regarding various sizes.
F_SIZE_OF_SMALL_STEPS = 1e-3
F_SIZE_OF_LARGE_STEPS = 1e-2
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

# Maximum # of steps when using Newton's method.
NEWTON_MAX_STEPS = 100