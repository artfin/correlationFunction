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

// hamiltonian input: theta(co2), phi, theta(h2), r
//#include "n2h2_hamiltonian_deriv.h"
// double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz);
//#include "co2h2_hamiltonian.h"

// dipole input: r, theta(co2), theta(h2), phi
//#include "n2h2_dipole_lr.h"
//double dipx(double R, double Theta1, double Theta2, double Phi);
//double dipy(double R, double Theta1, double Theta2, double Phi);
//double dipz(double R, double Theta1, double Theta2, double Phi);
//#include "n2h2_dip_derivatives.h"
//#include "n2h2_pot.h"
//#include "const_val.h"
//#include "constant.h"

//#include <Eigen/Dense>

//using namespace Eigen;

// potential: xr -- in angstroms
// pararead -- intializing procedure
// output of co2h2pes -- in cm^-1
// input: r, theta(co2), theta(h2), phi
extern "C" void __attribute__((stdcall)) pararead_(void);
extern "C" double co2h2pes_(double* rr, double* theta1, double* theta2, double* phi);


const int dimension = 4;

const double R0_num = 3.0;
const double R0_den = 4.0;

const double Rmax_num = 25.0;
const double Rmax_den = 25.0;

const double mu = 6632;
//const double mu = 38000;

const double a0 = 0.52917721067E-10; //in m
//double VOLUME = pow(Rmax_num,3.0)*4.0/3.0*M_PI;
const double VOLUME = 4.0/3.0*M_PI*pow(Rmax_den*a0,3.0)  ;


using namespace std;
using namespace std::placeholders;



int delta(int i, int j)
{
  return (i == j) ? 1 : 0;
}



double morse(double R)
{
	double val = sqrt(200)*(1- exp(-0.5*(R - 6.5)));
	return val*val - 200.0;
}




double density_full(double pR, double pTheta, double R, double Theta)
{
	// if ((x < 0.0) || (x > 10.0))
	// {
	// 	return x;
	// }
	// else
	// 	return 0.0;
	double H = pR*pR/2.0/mu+pTheta*pTheta/2.0/mu/R/R + ar_he_pot(R);

	//return exp(-H/3.166812E-6/300.0);
	//if ((H > 0) && (ar_he_pot(R) < 0))
	if ((H > 0))
	{
		return exp(-H/3.166812E-6/295.0);
	}
	else
		return 0.0;
}
 
void 
syst_2D (REAL t ,REAL *y, REAL *f)
{
   (void)(t); // avoid unused parameter warning 
  //double * par = (double *)params;

  

  double * out = new double[4];
  rhs_2D(out, y[0],y[1],y[2],y[3] );// R, pR, psi, ppsi
  


  f[0] = out[0];//dR/dt  
 // f[1]=f1(y,par)-(dpotdR);//d(pR)/dt 
  f[1]=out[1];
 
  f[2] = out[2];//d(psi)/dt 
  f[3] = out[3];//d(ppsi)/dt
 
  // f[4] = out[4];//d(phi)/dt
  // f[5] = out[5];//d(pPhi)/dt



  delete [] out;
  //delete [] par;

  
}


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
  fmax = 1000000; 


 
    if (time > 10.0)
    {
    	    fehler = gear4(&t0, xend, N, syst_2D, y0, epsabs,
                 epsrel, &h, fmax, &aufrufe);
    
    
	     if (fehler != 0) {
	     printf(" Gear4: error n° %d\n", 10 + fehler);
	   //  return 0;
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


	// double XJ, YJ, ZJ, Jac;

	// XJ = M_PI* (1 + xe*xe);
	// Jac = XJ;
	
	double R,   RJ;

	

	//double rnuc = sqrt(xe*xe);
	
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
	MPI_Init( &argc, &argv ); // ??

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	 

	double LTEMP = 300.0;
	double HTEMP = 300.0;
	double STEP = 5.0;

	// int number_of_temp = 36;
	// double Tem[number_of_temp] = {50,55,60,65,70,80,90,100,120,140,
	// 							150,160,170,180,200,220,240,250,260,270,
	// 							280,300,320,330,340,350,360,370,380,400,
	// 							420,440,450,460,480,500};

	 int number_of_temp = 1;
	 double Tem[number_of_temp] = {200};
	//int number_of_temp = 5;
	//double Tem[number_of_temp] = {380,400,420,450,480};

	

	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

	ofstream file;
	file.open( "second_moment_points.txt" );
	
	//hep::vegas_callback<double>(stop_after_precision(0.01));
	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	double TEMP;
	double TIME;
	double MAXTIME = 20000.0;
	double timstep = 400.0;

	//for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP )
	//for (int i = 0; i < number_of_temp; ++i)
	for (TIME = 8400.0; TIME <= MAXTIME; TIME += timstep )
	{
		//TIME = 0.0;
		
	 
		if ( rank == 0 )
		{	
			cout << "------------ correlation calculatin for time = "<< TIME << " ----------" << endl;
			//cout << "LTEMP = " << LTEMP << "; HTEMP = " << HTEMP << "; STEP = " << STEP << endl;
		}
		chrono::high_resolution_clock::time_point cycleStartTime = chrono::high_resolution_clock::now();
		
		auto numerator_integrand = bind( integrand_full, _1, TIME, true);
		auto denominator_integrand = bind( integrand_full, _1, TIME , false);

		auto numerator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(numerator_integrand, dimension),
			std::vector<std::size_t>(18, 6e4)
		);

		auto denominator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(denominator_integrand, dimension),
			std::vector<std::size_t>(12, 1e4)
		);

		auto numerator_result = hep::cumulative_result0( numerator_results.begin() + 2, numerator_results.end() );
		auto denominator_result = hep::cumulative_result0( denominator_results.begin() + 2, denominator_results.end() );

		double numerator_chi_square_dof = hep::chi_square_dof0( numerator_results.begin() + 2, numerator_results.end() );
		double denominator_chi_square_dof = hep::chi_square_dof0( denominator_results.begin() + 2, denominator_results.end() );

		double result = numerator_result.value()  / denominator_result.value()  * VOLUME ;
		//double result = numerator_result.value()*0.003614796373 ;

		if ( rank == 0 )
		{
			// cout << ">> Numerator; cumulative result: " << endl;
			//  cout << ">> N = " << numerator_result.calls() << "; I = " << numerator_result.value() << " +- " << numerator_result.error() << endl;
			//  cout << ">> chi^2/dof = " << numerator_chi_square_dof << endl << endl;

			// cout << ">> Denumerator; cumulative_result: " << endl;
			// cout << ">> N = " << denumerator_result.calls() << "; I = " << denumerator_result.value() << " +- " << denumerator_result.error() << endl;
			// cout << ">> chi^2/dof = " << denumerator_chi_square_dof << endl << endl;

			// file << setprecision(5) << Tem[i] << " " << setprecision(10) << numerator_result.value() << " " << denumerator_result.value() << " " << setprecision(10) << zero_moment << endl;
			//file <<  TIME*2.418884E-17 << " " <<  setprecision(10) << result << endl;
			file <<  TIME << " " <<  setprecision(10) << result << endl;
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