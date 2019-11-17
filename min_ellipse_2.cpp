#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/generators.h>

#include <CGAL/Min_ellipse_2.h>
#include "Min_ellipse_2_rublev.h"
#include "Min_ellipse_2_rublev_lp_solver.h"

#include <CGAL/Min_ellipse_2_traits_2.h>

#include <vector>
#include <cmath>
#include <ctime>
#include <cassert>
#include <string>

typedef  CGAL::Gmpq                       NT;
typedef  CGAL::Cartesian<NT>              K;
typedef  CGAL::Cartesian_d<NT>			  Kernel_d;
typedef  CGAL::Point_2<K>                 Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>  Traits;
typedef  CGAL::Point_d<Kernel_d>		  Point_d;

void
test_mve()
{
	int N, distribution, i;
	std::vector<long> sum, sum2;
	std::vector<double> M, M2;
	long prev_verbose_time, t;
	int iterations;
	std::list< Point> points;
	bool test_validness;
	std::string test_validness_s;

	const int num_of_methods = 4;

	sum.resize(num_of_methods, 0);
	sum2.resize(num_of_methods, 0);
	M.resize(num_of_methods, 0);
	M2.resize(num_of_methods, 0);

	std::cout << "Please, choose a distribution:\n" <<
		"\t1. in_square\n" <<
		"\t2. in_disk\n" <<
		"\t3. on_circle\n" <<
		"\t4. on_square\n";
	std::cin >> distribution;
	assert(distribution > 0 && distribution < 5);

	std::cout << "Please, enter N: ";
	std::cin >> N;
	assert(N > 0 && N < 1000000);

	std::cout << "Test validness ? (Y/n)";
	std::cin >> test_validness_s;
	if (test_validness_s.length() > 0 && (test_validness_s[0] == 'Y' || test_validness_s[0] == 'y'))
		test_validness = true;
	else
		test_validness = false;


//	std::cout << "        |             CGAL              |         rublev_fast\n";
//	std::cout << " Iter's |     median    |    std dev    |     median    |    std dev\n";
	std::cout << " Iter's |     rublev    |rublev_lp_solve|      CGAL     |rublev_convex_hull\n";

	iterations = 0;
	prev_verbose_time = 0;
	while ( true) {
		iterations ++;
		points.clear();

		if ( distribution == 1) {
			CGAL::Random_points_in_square_2< Point> generator(100.0);
			for ( i = 0; i < N; ++i)
				points.push_back( *( ++generator));
		}
		else if ( distribution == 2) {
			CGAL::Random_points_in_disc_2< Point> generator(100.0);
			for ( i = 0; i < N; ++i)
				points.push_back( *( ++generator));
		}
		else if ( distribution == 3) {
			CGAL::Random_points_on_circle_2< Point> generator(100.0);
			for ( i = 0; i < N; ++i)
				points.push_back( *( ++generator));
		}
		else if ( distribution == 4) {
			CGAL::Random_points_on_square_2< Point> generator(100.0);
			for ( i = 0; i < N; ++i)
				points.push_back( *( ++generator));
		}

		{
			t = clock();
			CGAL::Min_ellipse_2_rublev< Traits>  me( points.begin(), points.end(), false);
			t = clock() - t;
			sum[ 0] += t;
			sum2[ 0] += t*t;
			if (test_validness)
			{
				if ( !me.is_valid( false))
					me.is_valid( true);
			}
		}
		
		{
			t = clock();
			CGAL::Min_ellipse_2_rublev_lp_solver< Traits>  me( points.begin(), points.end());
			t = clock() - t;
			sum[ 1] += t;
			sum2[ 1] += t*t;
			if (test_validness)
			{
				if ( !me.is_valid( false))
					me.is_valid( true);
			}
		}

		{
			t = clock();
			CGAL::Min_ellipse_2< Traits>  me( points.begin(), points.end());
			t = clock() - t;
			sum[ 2] += t;
			sum2[ 2] += t*t;
			if (test_validness)
			{
				if ( !me.is_valid( false))
					me.is_valid( true);
			}
		}

		{
			t = clock();
			CGAL::Min_ellipse_2_rublev< Traits>  me( points.begin(), points.end(), true);
			t = clock() - t;
			sum[ 3] += t;
			sum2[ 3] += t*t;
			if (test_validness)
			{
				if ( !me.is_valid( false))
					me.is_valid( true);
			}
		}

		if (clock() - prev_verbose_time > 3000) {
			for (i = 0; i < num_of_methods; ++i)
			{
				M[ i] = ((double)sum[ i]/iterations) + 0.001;
				M2[ i] = sqrt( ((double) sum2[ i] / iterations) - M[ i] * M[ i]);
			}
				
			std::cout << std::setprecision(6) << " " << iterations;
				//<< " " << iterations << "\t| " << M[ 0] << "\t| " << M2[ 0] << "\t| " << M[ 1] << "\t| " << M2[ 1] << "\n";
			for (i = 0; i < num_of_methods; ++i)
				std::cout << std::setprecision(8) << "\t| " << M[ i];
			std::cout << "\n";
			std::cout.flush();
			prev_verbose_time = clock();
		}
	}
    return;
}

int
main( int, char**)
{
	test_mve();
}