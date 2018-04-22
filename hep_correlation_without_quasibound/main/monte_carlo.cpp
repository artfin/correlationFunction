#include <iostream>
#include <chrono>

#include <math.h>
#include <cstddef>
#include <vector>

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

//#include "hep/mc.hpp"
#include <hep/mc-mpi.hpp>
#include <gsl/gsl_sf_gamma.h>
#include "ar_he_dip_buryak_fit.h"
#include "ar_he_pes.h"
#include "diatomic_equations_2D.h"

/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "../jean-pierre/basis.h"         /*  for umleiten, fprintf, stderr, scanf,  */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include "../jean-pierre/vmblock.h"       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include "../jean-pierre/gear.h"          /*  for  gear4, gear_fehlertext            */
//#include "../../jean-pierre/backup/t_dgls.h"        /*  for  bsptyp, dgls_waehlen              */
/* ------------------------------------------------------------------ */

const int dimension = 4;

const double Rmin = 3.0;
const double Rmax = 20.0;

const double mu = 6632.09;

const double a0 = 0.52917721067E-10; //in m
const double VOLUME = 4.0/3.0*M_PI*pow(Rmax * a0, 3.0)  ;

class SplinedData
{
private:
    static std::vector<double> j_values;
    static std::vector<double> xmax_values;

public:
    SplinedData( std::string const & filename )
    {
        std::ifstream spline_file(filename);

        if ( !spline_file )
            throw std::invalid_argument("Can't find splined file.");

        const int MAX_LINE = 100;
        char line[MAX_LINE];

        while ( spline_file.getline(line, MAX_LINE) )
        {
            std::stringstream ss;
            ss.str(line);

            double j_value;
            ss >> j_value;

            double xmax_value;
            ss >> xmax_value;

            j_values.push_back( j_value );
            xmax_values.push_back( xmax_value );
        }

        spline_file.close();
    }

    double xmax( double j ) const
    {
        // нет максимума (либо прыгнули выше перегиба эффективного потенциала или слишком маленький угловой момент
        if  ( j < j_values[0] || j > j_values.back() )
        {
            //std::cout << "j less than" << j_values[0] << " or j more than " << j_values.back() << std::endl;
            return -1;
        }

        // нижняя граница для полученного значения j в сортированном массиве j_values
        std::vector<double>::iterator up_it = std::upper_bound( j_values.begin(), j_values.end(), j );
        size_t up = up_it - j_values.begin();

        // определяем верхнюю границу для полученного значения j в массиве j_values
        size_t low = up - 1;

        //std::cout << "input: " << j << "; low value: " << j_values[low] << "; low xmax: " << xmax_values[low] << std::endl;
        //std::cout << "input: " << j << "; up value: " << j_values[up] << "; up xmax: " << xmax_values[up] << std::endl;

        // если значение в списке
        if ( j == j_values[low] )
            return xmax_values[low];

        // иначе выдаем линейную аппроксимацию между найденными границами
        double k = (xmax_values[up] - xmax_values[low]) / (j_values[up] - j_values[low]);
        //std::cout << "k: " << k << std::endl;

        return xmax_values[up] - k * (j_values[up] - j);
    }
};
// инициализация статических векторов
std::vector<double> SplinedData::j_values;
std::vector<double> SplinedData::xmax_values;

// отбираем только пролетные траектории
double density_full(double pR, double pTheta, double R, double Theta, SplinedData const & sd )
{
    double H = pR * pR / 2.0 / mu + pTheta * pTheta / 2.0 / mu / R / R + ar_he_pot(R);

    // xmax = -1, если при данном J нет xmax
    // квазисвязанные (фактически связанные для двух атомов): R < Xmax(J)
    if ( R < sd.xmax(pTheta) )
        return 0.0;

    if ( H > 0 )
	{
        return std::exp(-H / 3.166812E-6 / 295.0);
	}
	else
		return 0.0;
}
 
void syst_2D (REAL t ,REAL *y, REAL *f)
{
   (void)(t); // avoid unused parameter warning 

  double * out = new double[4];
  rhs_2D(out, y[0],y[1],y[2],y[3] );// R, pR, psi, ppsi
  
  f[0] = out[0]; //dR/dt
  f[1] = out[1];
  f[2] = out[2]; //d(psi)/dt
  f[3] = out[3]; //d(ppsi)/dt
 
  delete [] out;
}

double correlation(double R, double pR, double psi, double ppsi, double time)
{
	REAL     epsabs = 1E-13;       /* absolute error bound                      */
	REAL     epsrel = 1E-13;       /* relative error bound                      */
  	REAL     t0 = 0.0;           /* left edge of integration interval         */
        
    REAL     h = 0.1;            /* initial, final step size                  */
    REAL     xend = time;         /* right edge of integration interval        */
    long     fmax = 1e8;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  	long     aufrufe;      /* actual number of function calls           */

	// 2d trajectories
  	int      N = 4;             //number of DEs in system                   
  	int      fehler;       /* error code from umleiten(), gear4()       */
  	int      i;            /* loop counter                              */
 
                        
    void     *vmblock = vminit();     /* List of dynamically allocated vectors     */
	/* out of memory? */ 
  	if (! vmcomplete(vmblock))  
	{       
    	printf("mgear: out of memory.\n");
    	throw std::runtime_error("out of memory!");
 	}
 
  	std::vector<REAL> y0 = {R, pR, psi, ppsi};

    double dipx0 = ar_he_dip_buryak_fit(R) * std::sin(psi);
    double dipy0 = ar_he_dip_buryak_fit(R) * std::cos(psi);
    double dipx, dipy;

    if (time > 0.0)
    {
        fehler = gear4(&t0, xend, N, syst_2D, &y0[0], epsabs, epsrel, &h, fmax, &aufrufe);

         if (fehler != 0)
         {
            printf(" Gear4: error n° %d\n", 10 + fehler);
            dipx = 0.0;
            dipy = 0.0;
         }
         else
         {
            dipx = ar_he_dip_buryak_fit(y0[0]) * std::sin(y0[2]);
            dipy = ar_he_dip_buryak_fit(y0[0]) * std::cos(y0[2]);
         }
    }
    else
    {
        dipx = dipx0;
        dipy = dipy0;
    }

    aufrufe = 0;
    vmfree( vmblock );

    return dipx0*dipx + dipy0*dipy;
}

double integrand_full( hep::mc_point<double> const& x, double time, SplinedData const & sd, bool is_num)
{	
	double pR = tan(M_PI * (x.point()[0] - 0.5));
	double pTheta =  tan(M_PI * (x.point()[1]/2.0));
	double Theta = x.point()[3] * M_PI * 2.0;
	double THJ = 2 * M_PI;
	double pRJ = M_PI* (1 + pR*pR);
	double pTJ = M_PI/2.0* (1 + pTheta*pTheta);
	double R  = x.point()[2] * (Rmax - Rmin) + Rmin;
    double RJ = Rmax - Rmin;

    double result;
    if (is_num)
        result = pTheta * correlation(R, pR, Theta, pTheta, time) * density_full(pR, pTheta, R, Theta, sd);
	else
        result = pTheta * density_full(pR, pTheta, R, Theta, sd);
	
    return result * RJ * THJ * pRJ * pTJ;
}

int main( int argc, char* argv[] )
{
    SplinedData sd("./j_xmax_values.txt");

	// initialize MPI
	MPI_Init( &argc, &argv ); 

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

    std::ofstream file( "correlation.txt" );
	
	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	double MINTIME = 0.0;
	double MAXTIME = 400.0;
	double timstep = 400.0;

	for (double TIME = MINTIME; TIME <= MAXTIME; TIME += timstep )
	{
		if ( rank == 0 )
            std::cout << "------------ correlation calculation for time = " << TIME << " ----------" << std::endl;
		
        std::chrono::high_resolution_clock::time_point cycleStartTime = std::chrono::high_resolution_clock::now();
		
        auto numerator_integrand = bind( integrand_full, std::placeholders::_1, TIME, sd, true);
        auto denominator_integrand = bind( integrand_full, std::placeholders::_1, TIME , sd, false);

		auto numerator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(numerator_integrand, dimension),
			std::vector<std::size_t>(10, 1e3)
		);

		auto denominator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(denominator_integrand, dimension),
			std::vector<std::size_t>(5, 1e5)
		);

		auto numerator_result = hep::cumulative_result0( numerator_results.begin() + 2, numerator_results.end() );
		auto denominator_result = hep::cumulative_result0( denominator_results.begin() + 2, denominator_results.end() );

		double numerator_chi_square_dof = hep::chi_square_dof0( numerator_results.begin() + 2, numerator_results.end() );
		double denominator_chi_square_dof = hep::chi_square_dof0( denominator_results.begin() + 2, denominator_results.end() );

		double result = numerator_result.value()  / denominator_result.value()  * VOLUME ;
		double error = numerator_result.error() / denominator_result.value() * VOLUME / result; // relative error of integration

		if ( rank == 0 )
		{
            file <<  TIME << " " <<  std::setprecision(10) << result << " " << error << std::endl;
            std::cout << "Final result: " << result << std::endl;
            std::chrono::high_resolution_clock::time_point cycleEndTime = std::chrono::high_resolution_clock::now();

            std::cout << "Cycle time: "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(cycleEndTime - cycleStartTime).count() / 1000.0
                      << " s" << std::endl;
		}
	}

	file.close();

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();

	if ( rank == 0 )
	{
        std::cout << std::endl << "Total time elapsed: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 1000.0
                  << " s" << std::endl;
	}
	
	MPI_Finalize();

	return 0;
}
