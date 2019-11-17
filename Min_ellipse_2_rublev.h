#ifndef CGAL_MIN_ELLIPSE_2_RUBLEV_H
#define CGAL_MIN_ELLIPSE_2_RUBLEV_H

#include <CGAL/Optimisation/basic.h>
#include <CGAL/Random.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <CGAL/ch_melkman.h>
#include <CGAL/Polygon_2.h>

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class Traits_ >
class Min_ellipse_2_rublev;

// Class interface
// ===============
template < class Traits_ >
class Min_ellipse_2_rublev {
  public:
    // types
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Ellipse                  Ellipse;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;
    
  private:
    // private data members
    Traits       tco;                           // traits class object
    std::list<Point>  points;                   // doubly linked list of points
    int          n_support_points;              // number of support points
    Point*       support_points;                // array of support points
	
	bool use_convex_hull_heuristic;
	int n_convex_hull_points;
	Point* convex_hull;    

    // copying and assignment not allowed!
    Min_ellipse_2_rublev( const Min_ellipse_2_rublev<Traits_>&);
    Min_ellipse_2_rublev<Traits_>& operator = ( const Min_ellipse_2_rublev<Traits_>&);

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
	void update_ellipse_old( Point add )
	{
		  int i;
		  Min_ellipse_2< Traits_ > algo( add);
		  algo.insert( support_points + 0, support_points + n_support_points);
		  n_support_points = algo.number_of_support_points();
		  for (i = 0; i < n_support_points; ++i)
			  support_points[i] = algo.support_point(i);
		  tco = algo.traits();
	} 

	inline void
		compute_ellipse( )
    {
        switch ( n_support_points) {
          case 5:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2],
                             support_points[ 3],
                             support_points[ 4]);
            break;
          case 4:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2],
                             support_points[ 3]);
            break;
          case 3:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2]);
            break;
          case 2:
            tco.ellipse.set( support_points[ 0], support_points[ 1]);
            break;
          case 1:
            tco.ellipse.set( support_points[ 0]);
            break;
          case 0:
            tco.ellipse.set( );
            break;
          default:
            CGAL_optimisation_assertion( ( n_support_points >= 0) &&
                                         ( n_support_points <= 5) ); }
    }

	inline void 
		me_recursive_support ( Point* first, const Point* last, int n_sp)
    {
        // compute ellipse through support points
        n_support_points = n_sp;
        compute_ellipse();
        if ( n_sp == 5) return;
    
        // test first n points
		Point* point_iter;
        for ( 
			point_iter = first; 
			last != point_iter; )
		{
            const Point& p = *point_iter;
    
            // p not in current ellipse?
            if ( has_on_unbounded_side( p))
			{
				// recursive call with p as additional support point
                support_points[ n_sp] = p;
                me_recursive_support( first, point_iter, n_sp+1);
			}
            else
                ++point_iter;
		}
    }

	inline void
		update_ellipse( Point add )
	{
		Point points[5];
		int n_sp = n_support_points;
		for (int i = 0; i < n_support_points; ++i)
			points[ i] = support_points[ i];
		
		support_points[ 0] = add;

		me_recursive_support( points + 0, points + n_sp, 1);
	} 

	inline void
		compute_convex_hull ()
	{
		std::vector< Point > basis;
		ch_melkman ( suppo	rt_points + 0, support_points + n_support_points , std::back_inserter( basis), Point::R ());

		n_convex_hull_points = basis.size();
		for (int i = 0; i < basis.size(); ++i)
			convex_hull[ i] = basis[ i];
	}

	inline void
    me_with_convex_hull_heuristic( )
    {
		std::list< Point > 
			iter_points[ 2];
		typename std::list< Point >::iterator last, point_iter;
		int cur_iter = 0;
		bool outlier = true;

		compute_ellipse();
		std::copy( points.begin(), points.end(), std::back_inserter( iter_points[ cur_iter]));

		compute_convex_hull();

		while ( outlier)
		{
			outlier = false;
			last = iter_points[ cur_iter].end();

			for ( point_iter = iter_points[ cur_iter].begin(); point_iter != last; point_iter++)
			{
				const Point& p = *point_iter;

				if ( bounded_side_2( convex_hull + 0, convex_hull + n_convex_hull_points, p, Point::R()) != ON_BOUNDED_SIDE)
				{
					iter_points[ 1-cur_iter].push_back( p);

					if ( tco.ellipse.has_on_unbounded_side( p))
					{
						update_ellipse( p);
						outlier = true;
						compute_convex_hull();
					}
				}
			}

			iter_points[ cur_iter].clear();
			cur_iter = 1-cur_iter;
		}
	}
	
		
	inline void
    me ( )
    {
		typename std::list< Point >::iterator last = points.end(), it;
		bool outlier = true;

		n_support_points = 0;
		tco.ellipse.set();

		while ( outlier)
		{
			outlier = false;
			for ( it = points.begin(); it != last; it++)
			{
				const Point& p = *it;
				if ( tco.ellipse.has_on_unbounded_side( p))
				{
					update_ellipse( p);
					outlier = true;
				}
			}
		}
    }

  public:
    // Constructors
    // ------------
    // STL-like constructor (member template)
    template < class InputIterator >
    Min_ellipse_2_rublev( InputIterator first,
                   InputIterator last,
				   bool use_convex_hull_heuristic,
                   const Traits& traits    = Traits())
            : tco( traits)
        {
			this->use_convex_hull_heuristic = use_convex_hull_heuristic;

            // allocate support points' array
            support_points = new Point[ 5];
			n_support_points = 0;

			convex_hull = new Point[ 5];
			n_convex_hull_points = 0;
    
            // range of points not empty?
            if ( first != last) {    
                // store points
                std::copy( first, last, std::back_inserter( points));
			}
    
            // compute me
			if ( use_convex_hull_heuristic)
			{
				me_with_convex_hull_heuristic( );
			}
			else
			{
				me( );
			}
        }

    // default constructor
    inline
    Min_ellipse_2_rublev( const Traits& traits = Traits())
        : tco( traits), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		n_support_points = 0;

		convex_hull = new Point[ 5];
		n_convex_hull_points = 0;
    
        // initialize ellipse
        tco.ellipse.set();
    
        CGAL_optimisation_postcondition( is_empty());
    }
    
    // constructor for one point
    inline
    Min_ellipse_2_rublev( const Point& p, const Traits& traits = Traits())
        : tco( traits), points( 1, p), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		convex_hull = new Point[ 5];
    
        // initialize ellipse
		n_support_points = 1;
        support_points[ 0] = p;
        tco.ellipse.set( p);

		n_convex_hull_points = 1;
		convex_hull[ 0] = p;
    
        CGAL_optimisation_postcondition( is_degenerate());
    }
    
    // constructor for two points
    // This was const Point& but then Intel 7.0/.net2003 messes it up 
    // with the constructor taking an iterator range
    inline
    Min_ellipse_2_rublev( Point p1, Point p2,
                   const Traits& traits = Traits())
        : tco( traits), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		n_support_points = 0;

		convex_hull = new Point[ 5];
		n_convex_hull_points = 0;
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
    
        // compute me
        me_with_convex_hull_heuristic( );
    
        CGAL_optimisation_postcondition( is_degenerate());
    }
    
    // constructor for three points
    inline
    Min_ellipse_2_rublev( const Point& p1, const Point& p2, const Point& p3,
                   const Traits& traits = Traits())
        : tco( traits), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		n_support_points = 0;

		convex_hull = new Point[ 5];
		n_convex_hull_points = 0;
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
    
        // compute me
        me_with_convex_hull_heuristic();
    }
    
    // constructor for four points
    inline
    Min_ellipse_2_rublev( const Point& p1, const Point& p2,
                   const Point& p3, const Point& p4,
                   const Traits& traits = Traits())
        : tco( traits), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		n_support_points = 0;

		convex_hull = new Point[ 5];
		n_convex_hull_points = 0;
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
        points.push_back( p4);
    
        // compute me
        me_with_convex_hull_heuristic();
    }
    
    // constructor for five points
    inline
    Min_ellipse_2_rublev( const Point& p1, const Point& p2, const Point& p3,
                   const Point& p4, const Point& p5,
                   const Traits& traits = Traits())
        : tco( traits), use_convex_hull_heuristic( true)
    {
        // allocate support points' array
        support_points = new Point[ 5];
		n_support_points = 0;

		convex_hull = new Point[ 5];
		n_convex_hull_points = 0;
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
        points.push_back( p4);
        points.push_back( p5);
    
        // compute me
        me_with_convex_hull_heuristic();
    }
    

    // Destructor
    // ----------
    inline
    ~Min_ellipse_2_rublev( )
    {
        // free support points' array
        delete[] support_points;
		delete[] convex_hull;
    }

    // Modifiers
    // ---------
    void
    insert( const Point& p)
    {
		if ( !use_convex_hull_heuristic)
		{
			use_convex_hull_heuristic = true;
			compute_convex_hull();
		}

		// p not in current support convex hull?
		if ( bounded_side_2( convex_hull + 0, convex_hull + n_convex_hull_points, p, Point::R()) != ON_BOUNDED_SIDE)
		{
			if ( has_on_unbounded_side( p))
			{
				points.push_front( p);
				support_points[0] = p;
				n_support_points = 1;

				me_with_convex_hull_heuristic();
			}
				// append p to the end of the list
				points.push_back( p);
		}
		else
			// append p to the end of the list
            points.push_back( p);
    }

    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
    {
        for ( ; first != last; ++first)
            insert( *first);
    }

    void
    clear( )
    {
        points.erase( points.begin(), points.end());
        n_support_points = 0;
		n_convex_hull_points = 0;

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

#endif // CGAL_MIN_ELLIPSE_2_RUBLEV_H

// ===== EOF =================================================================