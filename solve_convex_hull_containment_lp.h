// example: function to check whether a point is in the convex 
// hull of other points; this version uses a maker
#include <boost/config.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Kernel_traits.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// unary function to get homogeneous begin-iterator of point
template <class Point_d>
struct Homogeneous_begin  {
  typedef typename Point_d::Homogeneous_const_iterator result_type;
  result_type operator() (const Point_d& p) const {
    return p.homogeneous_begin();
  }
};

// function to test whether point is in the convex hull of other points;
// the type ET is an exact type used for the computations
template <class Point_d, class RandomAccessIterator, class ET>
CGAL::Quadratic_program_solution<ET>
solve_convex_hull_containment_lp (const Point_d& p,
				  RandomAccessIterator begin,
				  RandomAccessIterator end, const ET& dummy)
{
  // construct program and solve it
  return CGAL::solve_nonnegative_linear_program
    (CGAL::make_nonnegative_linear_program_from_iterators
     (end-begin,                                                           // n
      p.dimension()+1,                                                     // m
      boost::transform_iterator
      <Homogeneous_begin<Point_d>, RandomAccessIterator>(begin),           // A
      typename Point_d::Homogeneous_const_iterator (p.homogeneous_begin()),// b
      CGAL::Const_oneset_iterator<CGAL::Comparison_result>(CGAL::EQUAL),   // ~
      CGAL::Const_oneset_iterator
      <typename CGAL::Kernel_traits<Point_d>::Kernel::RT> (0)), dummy);    // c
} 

// I've been puzzled with the example provided in CGAL for testing whether the point is inside the convex hull of point set.
// The main issue was that Point_2 (in comparison with Point_d) doesn't contain an Homogeneous_const_iterator type. 
// I tried to use something like Cartesian_const_iterator instead, but then some runtime assertion is raised. 
// So I've came up with dumb approach simply converting all available data to what is being used in example.
// I also have hard times every time I need to debug something. So the process is not so fast as I wish.
template <class Point, class RandomAccessIterator, class Point_d, class ET> 
bool is_in_convex_hull (const Point& p,
			RandomAccessIterator begin,
			RandomAccessIterator end,
			const Point_d& dummy,
			const ET& dummy2)
{
	Point_d point = Point_d(p.x(), p.y());
	std::vector<Point_d> points;
	for (std::vector<Point>::iterator it = begin; it != end; ++it)
		points.push_back(Point_d(it->x(), it->y()));

	CGAL::Quadratic_program_solution<ET> s =
		solve_convex_hull_containment_lp (point, points.begin(), points.end(), ET(0));

    return !s.is_infeasible();
}
