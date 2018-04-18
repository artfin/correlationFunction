#include <iostream>
#include <chrono>

#include <math.h>
#include <cstddef>
#include <vector>

#include <fstream>
#include <iomanip>
#include <stdexcept>

//#include "hep/mc.hpp"
#include <hep/mc-mpi.hpp>
#include <gsl/gsl_sf_gamma.h>
#include "ar_he_dip_buryak_fit.h"
#include "ar_he_pes.h"
#include "diatomic_equations_2D.h"


//#include "const_values.h"
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

const double R0_num = 3.0;
const double R0_den = 4.0;

const double Rmax_num = 25.0;
const double Rmax_den = 25.0;

const double mu = 6632;
//const double mu = 38000;

const double a0 = 0.52917721067E-10; //in m
const double VOLUME = 4.0/3.0*M_PI*pow(Rmax_den*a0,3.0)  ;


using namespace std;
using namespace std::placeholders;

double density_full(double pR, double pTheta, double R, double Theta)
{
	double H = pR*pR/2.0/mu+pTheta*pTheta/2.0/mu/R/R + ar_he_pot(R);

	if ((H > 0))
	{
		return exp(-H/3.166812E-6/295.0);
	}
	else
		return 0.0;
}
 
void syst_2D (REAL t ,REAL *y, REAL *f)
{
   (void)(t); // avoid unused parameter warning 

  double * out = new double[4];
  rhs_2D(out, y[0],y[1],y[2],y[3] );// R, pR, psi, ppsi
  
  f[0] = out[0];//dR/dt  
 // f[1]=f1(y,par)-(dpotdR);//d(pR)/dt 
  f[1]=out[1];
 
  f[2] = out[2];//d(psi)/dt 
  f[3] = out[3];//d(ppsi)/dt
 
  delete [] out;
}

struct stop_after_precision
{
	stop_after_precision( double abs_error, double rel_error )
		: abs_error(abs_error), rel_error(rel_error)
	{
	}

	double abs_error;
	double rel_error;

	bool operator()(std::vector<hep::vegas_result<double>> const& r )
	{
		// print the results obtained so far
		hep::vegas_verbose_callback<double>(r);

		// compute cumulative result...
		auto const result = hep::accumulate<hep::weighted_with_variance>(r.begin(), r.end() );

		// check absolute and relative errors
		if ( result.error() < abs_error )
		{
			std::cout << " absolute error " << result.error() << " is smaller than the limit " << abs_error << std::endl;

			return false;
		}

		if ( result.error() < rel_error * result.value() )
		{
			std::cout << " relative error " << (result.error() / result.value()) << " is smaller than the limit " << rel_error << std::endl;
			
			return false;
		}

		return true;
	}
};


double correlation(double R, double pR, double psi, double ppsi, double time)
{

  REAL     epsabs;       /* absolute error bound                      */
  REAL     epsrel;       /* relative error bound                      */
  REAL     t0;           /* left edge of integration interval         */
  REAL     *y0;          /* [0..n-1]-vector: initial value, approxim. */
        
  REAL     h;            /* initial, final step size                  */
  REAL     xend;         /* right edge of integration interval        */
  long     fmax;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  long     aufrufe;      /* actual number of function calls           */

  int      N;             //number of DEs in system                   
  int      fehler;       /* error code from umleiten(), gear4()       */
  int      i;            /* loop counter                              */
 
                        
  void     *vmblock;     /* List of dynamically allocated vectors     */

 // FILE * trajectory_file = fopen("correlation_tests.txt","w");

/* -------------------- read input  -------------------- */

  //N = 6; //3D
  N = 4;//2D
 
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
 
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    throw std::runtime_error("out of memory!");
   // return 0;
  }
 
  //double end = time;
   

  double dipx, dipy, dipz;
  double dipx0, dipy0, dipz0;

 
  epsabs = 1E-13;
  epsrel = 1E-13;

  t0 = 0.0;
 
  y0[0] = R;
  y0[1] = pR;
  y0[2] = psi;
  y0[3] = ppsi;


  dipx0 = ar_he_dip_buryak_fit(R)*sin(psi);
  dipy0 = ar_he_dip_buryak_fit(R)*cos(psi);
   
  h =0.1;
  xend = time;
  fmax = 1e6; 
    
  if (time > 10.0)
  {
     fehler = gear4(&t0, xend, N, syst_2D, y0, epsabs, epsrel, &h, fmax, &aufrufe);
	 
	 if (fehler != 0) {
	     printf(" Gear4: error n° %d\n", 10 + fehler);
	     dipx = 0.0;
	     dipy = 0.0;
	}
	else
	{
	  	dipx = ar_he_dip_buryak_fit(y0[0])*sin(y0[2]);
	 	dipy = ar_he_dip_buryak_fit(y0[0])*cos(y0[2]);
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

double integrand_full( hep::mc_point<double> const& x, double time, bool is_num)
{	
	double pR = tan(M_PI * (x.point()[0] - 0.5));
	//double pTheta = tan(M_PI * (x.point()[1] - 0.5));
	double pTheta =  tan(M_PI * (x.point()[1]/2.0));

	//double R  = x.point()[2] * (Rmax_num-R0_num) + R0_num;
	double Theta = x.point()[3] * M_PI * 2.0;
	

	//double RJ = (Rmax_num-R0_num);
	double THJ = 2 * M_PI;
	double pRJ = M_PI* (1 + pR*pR);
	//double pTJ = M_PI* (1 + pTheta*pTheta);
	double pTJ = M_PI/2.0* (1 + pTheta*pTheta);

	double R,   RJ;

	double result; 
	if (is_num)
	{
		R  = x.point()[2] * (Rmax_num-R0_num) + R0_num;
		//result = pTheta* ar_he_dip_buryak_fit(R)*ar_he_dip_buryak_fit(R)*density_full(pR, pTheta, R, Theta);
		result = pTheta* correlation(R, pR, Theta, pTheta, time)*density_full(pR, pTheta, R, Theta);
		//result = density_full(pR, pTheta, R, Theta);
		RJ = (Rmax_num-R0_num);
	}
	else
	{
		R  = x.point()[2] * (Rmax_den-R0_den) + R0_den;
		result = pTheta* density_full(pR, pTheta, R, Theta);
		RJ = (Rmax_den-R0_den);
	}
	return result*RJ*THJ*pRJ*pTJ;
}

int main( int argc, char* argv[] )
{
	// initialize MPI
	MPI_Init( &argc, &argv ); 

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

	ofstream file;
	file.open( "correlation.txt" );
	
	// понятия не имею какая должна быть абсолютная ошибка, будет считаться пока относительная не станет < 1e-3 = 0.1% 
	// надо выставить циклов с запасом
	hep::vegas_callback<double>( stop_after_precision(1e-50, 1e-3) ); 
	

	double MINTIME = 0.0;
	double MAXTIME = 2000.0;
	double timstep = 400.0;

	for (double TIME = MINTIME; TIME <= MAXTIME; TIME += timstep )
	{
		if ( rank == 0 )
			cout << "------------ correlation calculation for time = " << TIME << " ----------" << endl;
		
		chrono::high_resolution_clock::time_point cycleStartTime = chrono::high_resolution_clock::now();
		
		auto numerator_integrand = bind( integrand_full, _1, TIME, true);
		auto denominator_integrand = bind( integrand_full, _1, TIME , false);

		auto numerator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(numerator_integrand, dimension),
			std::vector<std::size_t>(15, 1e4)
		);

		auto denominator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(denominator_integrand, dimension),
			std::vector<std::size_t>(20, 1e5)
		);

		auto numerator_result = hep::accumulate<hep::weighted_with_variance>(
        	numerator_results.begin() + 1, numerator_results.end() );
	
		auto denominator_result = hep::accumulate<hep::weighted_with_variance>(
			denominator_results.begin() + 1, denominator_results.end() );

		double result = numerator_result.value()  / denominator_result.value()  * VOLUME ;
		double error = numerator_result.error() / denominator_result.value() * VOLUME / result; // relative error of integration

		if ( rank == 0 )
		{
			file <<  TIME << " " <<  setprecision(10) << result << " " << error << endl;
			cout << "Final result: " << result<< endl;
			chrono::high_resolution_clock::time_point cycleEndTime = chrono::high_resolution_clock::now();

			cout << "Cycle time: " << chrono::duration_cast<chrono::milliseconds>(cycleEndTime - cycleStartTime).count() / 1000.0 << " s" << endl; 
		}
	}

	file.close();

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

	if ( rank == 0 )
	{
		cout << endl << "Total time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 << " s" << endl;
	}
	
	MPI_Finalize();

	return 0;
}


