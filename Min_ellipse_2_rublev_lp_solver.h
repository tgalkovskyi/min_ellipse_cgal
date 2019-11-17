#ifndef CGAL_MIN_ELLIPSE_2_RUBLEV_LP_SOLVER_H
#define CGAL_MIN_ELLIPSE_2_RUBLEV_LP_SOLVER_H

#include <CGAL/Optimisation/basic.h>
#include <CGAL/Random.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian_d.h>

#include "solve_convex_hull_containment_lp.h"

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class Traits_ >
class Min_ellipse_2_rublev_lp_solver;

// Class interface
// ===============
template < class Traits_ >
class Min_ellipse_2_rublev_lp_solver {
  public:
    // types
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Ellipse                  Ellipse;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;
	typedef typename  Traits_::K::FT NT;
	typedef typename  Cartesian_d<NT>::Point_d Point_d;
    
  private:
    // private data members
    Traits       tco;                           // traits class object
    std::list<Point>  points;                   // doubly linked list of points
    int          n_support_points;              // number of support points
    Point*       support_points;                // array of support points
    

    // copying and assignment not allowed!
    Min_ellipse_2_rublev_lp_solver( const Min_ellipse_2_rublev_lp_solver<Traits_>&);
    Min_ellipse_2_rublev_lp_solver<Traits_>& operator = ( const Min_ellipse_2_rublev_lp_solver<Traits_>&);

// ============================================================================

// Class implementation
// ====================

  public:
    // Access functions and predicates
    // -------------------------------
    // #points and #support points
    inline
    int
    number_of_points( ) const
    {
        return( points.size());
    }
    
    inline
    int
    number_of_support_points( ) const
    {
        return( n_support_points);
    }

    // is_... predicates
    inline
    bool
    is_empty( ) const
    {
        return( number_of_support_points() == 0);
    }
    
    inline
    bool
    is_degenerate( ) const
    {
        return( number_of_support_points() <  3);
    }

    // access to points and support points
    inline
    Point_iterator
    points_begin( ) const
    {
        return( points.begin());
    }
    
    inline
    Point_iterator
    points_end( ) const
    {
        return( points.end());
    }
    
    inline
    Support_point_iterator
    support_points_begin( ) const
    {
        return( support_points);
    }
    
    inline
    Support_point_iterator
    support_points_end( ) const
    {
        return( support_points+n_support_points);
    }
    
    // random access for support points
    inline
    const Point&
    support_point( int i) const
    {
        CGAL_optimisation_precondition( (i >= 0) &&
                                        (i <  number_of_support_points()));
        return( support_points[ i]);
    }
    // ellipse
    inline
    const Ellipse&
    ellipse( ) const
    {
        return( tco.ellipse);
    }
    
    // in-ellipse test predicates
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( tco.ellipse.bounded_side( p));
    }
    
    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_bounded_side( p));
    }
    
    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( tco.ellipse.has_on_boundary( p));
    }
    
    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_unbounded_side( p));
    }

  private:
    // Private member functions
    // ------------------------
	inline
	void compute_ellipse( Point add )
	  {
		  int i;
		  Min_ellipse_2< Traits_ > algo( add);
		  algo.insert( support_points + 0, support_points + n_support_points);
		  n_support_points = algo.number_of_support_points();
		  for (i = 0; i < n_support_points; ++i)
			  support_points[i] = algo.support_point(i);
		  tco = algo.traits();
	  } 

    void
    me( )
    {
		std::list< Point > 
			iter_points[ 2], 
			convex_hull;
		typename std::list< Point >::iterator last;
		int cur_iter = 0;
		bool outlier = true;

		n_support_points = 0;
		tco.ellipse.set();
		std::copy( points.begin(), points.end(), std::back_inserter( iter_points[ cur_iter]));

		while ( outlier) {
			outlier = false;
			last = iter_points[ cur_iter].end();

			for ( typename std::list<Point>::iterator point_iter = iter_points[ cur_iter].begin(); point_iter != last; point_iter++) {
				const Point& p = *point_iter;

				if ( !is_in_convex_hull( p, support_points + 0, support_points + n_support_points, Point_d(), NT(0)))
				{
					iter_points[ 1-cur_iter].push_back( p);

					if ( tco.ellipse.has_on_unbounded_side( p)) {
						compute_ellipse( p);
						outlier = true;
					}
				}
			}

			iter_points[ cur_iter].clear();
			cur_iter = 1-cur_iter;
		}
    }

  public:
    // Constructors
    // ------------
    // STL-like constructor (member template)
    template < class InputIterator >
    Min_ellipse_2_rublev_lp_solver( InputIterator first,
                   InputIterator last,
                   const Traits& traits    = Traits())
            : tco( traits)
        {
            // allocate support points' array
            support_points = new Point[ 5];
    
            // range of points not empty?
            if ( first != last) {    
                // store points
                std::copy( first, last, std::back_inserter( points));
			}
    
            // compute me
			me();
        }

    // Destructor
    // ----------
    inline
    ~Min_ellipse_2_rublev_lp_solver( )
    {
        // free support points' array
        delete[] support_points;
    }

    void
    clear( )
    {
        points.erase( points.begin(), points.end());
        n_support_points = 0;
        tco.ellipse.set();
    }
    

    // Validity check
    // --------------
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        using namespace std;
    
        CGAL::Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL::Min_ellipse_2_rublev<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;
    
        // containment check (a)
        verr << "  a) containment check..." << flush;
        Point_iterator point_iter;
        for ( point_iter  = points_begin();
              point_iter != points_end();
              ++point_iter)
            if ( has_on_unbounded_side( *point_iter))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "ellipse does not contain all points"));
        verr << "passed." << endl;
    
        // support set checks (b)+(c) (not yet implemented)
        
        // alternative support set check
        verr << "  +) support set check..." << flush;
        Support_point_iterator support_point_iter;
        for ( support_point_iter  = support_points_begin();
              support_point_iter != support_points_end();
              ++support_point_iter)
            if ( ! has_on_boundary( *support_point_iter))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "ellipse does not have all \
                             support points on the boundary"));
        verr << "passed." << endl;
    
        verr << "  object is valid!" << endl;
        return( true);
    }

    // Miscellaneous
    // -------------
    inline
    const Traits&
    traits( ) const
    {
        return( tco);
    }
};

//TODO: IO functions needs to be added

CGAL_END_NAMESPACE

#endif // CGAL_MIN_ELLIPSE_2_RUBLEV_LP_SOLVER_H

// ===== EOF =================================================================