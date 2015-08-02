// Computing intersection points among curves using the sweep line.
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <vector>
#include <list>

#include "cgal_intersection.h"

typedef CGAL::Cartesian<double>                         Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Polyline_traits;
typedef Polyline_traits::Point_2                        Point;
typedef Polyline_traits::Curve_2                        Polyline;

extern "C" int find_intersections_of_curves(
    coordinate *curve_1, long curve_1_size,
    coordinate *curve_2, long curve_2_size,
    coordinate *intersections, int max_num_of_intersections)
{
    Polyline_traits polyline_traits;
    Polyline_traits::Construct_curve_2 construct_polyline =
        polyline_traits.construct_curve_2_object();

    int num_of_curves = 2;
    coordinate *curves[] = {curve_1, curve_2};
    long size_of_curve[] = {curve_1_size, curve_2_size};
    std::vector<Polyline> polylines;

    for(int n = 0; n < num_of_curves; n++){
        std::vector<Point> points;
        Point p_0 = Point(curves[n][0].x, curves[n][0].y);
        points.push_back(p_0);
        for(long i = 1; i < size_of_curve[n]; i++){
            Point p_i = Point(curves[n][i].x, curves[n][i].y);
            if(p_i != points.back()){
                points.push_back(p_i);
            }
        }
        polylines.push_back(construct_polyline(points.begin(), points.end()));
    }

    std::list<Point> intersection_points;
    CGAL::compute_intersection_points(polylines.begin(), polylines.end(),
                                      std::back_inserter(intersection_points),
                                      false, polyline_traits);

    int num_of_intersections = intersection_points.size();
    Point intersection_point;
    for (int j = 0; j < num_of_intersections; j++){
        if (j == max_num_of_intersections)
            break;
        intersection_point = intersection_points.front();
        intersections[j].x = intersection_point.x();
        intersections[j].y = intersection_point.y();
        intersection_points.pop_front();
    }

    return num_of_intersections;
}
