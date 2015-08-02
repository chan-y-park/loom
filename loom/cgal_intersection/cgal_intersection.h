typedef struct {double x; double y;} coordinate;

int cgal();

extern "C" int find_intersections_of_curves(
    coordinate *curve_1, long curve_1_size,
    coordinate *curve_2, long curve_2_size,
    coordinate *intersections, int max_num_of_intersections);
