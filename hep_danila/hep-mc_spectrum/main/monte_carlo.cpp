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


double correlation(double R, double pR, double psi, double ppsi, double time)
{
	REAL     epsabs = 1E-13;       /* absolute error bound                      */
	REAL     epsrel = 1E-13;       /* relative error bound                      */
  	REAL     t0 = 0.0;           /* left edge of integration interval         */
        
  	REAL     h;            /* initial, final step size                  */
  	REAL     xend;         /* right edge of integration interval        */
  	long     fmax;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  	long     aufrufe;      /* actual number of function calls           */

	// 2d trajectories
  	int      N = 4;             //number of DEs in system                   
  	int      fehler;       /* error code from umleiten(), gear4()       */
  	int      i;            /* loop counter                              */
 
                        
  	void     *vmblock;     /* List of dynamically allocated vectors     */

  	vmblock = vminit();                 /* initialize storage */
	
	/* out of memory? */ 
  	if (! vmcomplete(vmblock))  
	{       
    	printf("mgear: out of memory.\n");
    	throw std::runtime_error("out of memory!");
 	}
 
  	std::vector<REAL> y0 = {R, pR, psi, ppsi};
	std::vector<REAL> y0_back{R, -pR, psi, -ppsi};

	std::vector<double> dipx_forward, dipx_backward;
	std::vector<double> dipy_forward, dipy_backward;

	dipx_forward.push_back( ar_he_dip_buryak_fit(R) * sin(psi) );
	dipx_backward.push_back( ar_he_dip_buryak_fit(R) * sin(psi) );
	
	dipy_forward.push_back( ar_he_dip_buryak_fit(R) * cos(psi) );
	dipy_backward.push_back( ar_he_dip_buryak_fit(R) * cos(psi) );

  	h = 0.1;
	xend = time;
  	fmax = 1e8; 

	if ( time == 0.0 )
		return ar_he_dip_buryak_fit(R) * ar_he_dip_buryak_fit(R);

   	const int MaxTrajectoryLength = 1000;	
	const double sampling_time = time;

  	for ( int step_counter = 0; y0[0] < Rmax; step_counter++, xend += time )  
  	{
   		fehler = gear4(&t0, xend, N, syst_2D, &y0[0], epsabs, epsrel, &h, fmax, &aufrufe);
	 
	 	if (fehler != 0) 
	 	{
	   		printf(" Gear4: error n° %d\n", 10 + fehler);
			return 0;
		}

		if ( step_counter == MaxTrajectoryLength )
			break;
	
		dipx_forward.push_back( ar_he_dip_buryak_fit(y0[0]) * sin(y0[2]) );
	 	dipy_forward.push_back( ar_he_dip_buryak_fit(y0[0]) * cos(y0[2]) );
 	}

	t0 = 0.0;
	xend = time;
	for ( int step_counter = 0; y0_back[0] < Rmax; step_counter++, xend += time )
	{
		fehler = gear4( &t0, xend, N, syst_2D, &y0_back[0], epsabs, epsrel, &h, fmax, &aufrufe );

		if ( fehler != 0 )
		{
			std::cout << "Gear4: error n " << 10 + fehler << std::endl;
			return 0;
		}

		if ( step_counter == MaxTrajectoryLength )
			break;

		dipx_backward.push_back( ar_he_dip_buryak_fit(y0_back[0]) * sin(y0_back[2]) );
		dipy_backward.push_back( ar_he_dip_buryak_fit(y0_back[0]) * sin(y0_back[2]) );
	}

	double corr = 0.0;
	if ( dipx_forward.size() > 1 )
		for ( size_t k = 1; k < dipx_forward.size(); ++k )
			corr += ( dipx_forward[k - 1] * dipx_forward[k] + dipy_forward[k - 1] * dipy_forward[k] );

	if ( dipx_backward.size() > 1 )
		for ( size_t k = 1; k < dipx_backward.size(); ++k )
			corr += ( dipx_backward[k - 1] * dipx_backward[k] + dipy_backward[k - 1] * dipy_backward[k] );

	corr /= ( dipx_forward.size() + dipx_backward.size() - 2 );

    aufrufe = 0;
 	vmfree( vmblock );
  
	return corr; 
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

	double result; 
	double R  = x.point()[2] * (Rmax - Rmin) + Rmin;
	double RJ = Rmax - Rmin; 
	if (is_num)
		result = pTheta * correlation(R, pR, Theta, pTheta, time) * density_full(pR, pTheta, R, Theta);
	
	else
		result = pTheta * density_full(pR, pTheta, R, Theta);
	
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
	
	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	double MINTIME = 0.0;
	double MAXTIME = 400.0;
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
