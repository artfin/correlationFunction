#include <iostream>
#include <chrono>

#include <math.h>
#include <cstddef>
#include <vector>

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include <fftw3.h>

#include "hep/mc-mpi.hpp"
#include <gsl/gsl_sf_gamma.h>
#include "ar_he_dip_buryak_fit.hpp"
#include "ar_he_pes.hpp"
#include "diatomic_equations_3d.hpp"
#include "constants.hpp"

/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "basis.hpp"         /*  for umleiten, fprintf, stderr, scanf,  */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include "vmblock.hpp"       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include "gear.hpp"          /*  for  gear4, gear_fehlertext            */
/* ------------------------------------------------------------------ */

const double Temperature = 295.0;

const double R0_min = 3.0;
const double R0_max = 25.0;

const double mu = 6632.039; // He-Ar

const double VOLUME = 4.0/3.0*M_PI*std::pow(R0_max * constants::ALU, 3.0);

const int MaxTrajectoryLength = 8092;
const int FFTMaxLength = MaxTrajectoryLength / 2 + 1;
double sampling_time = 100; // a.t.u.

double kT = constants::BOLTZCONST * Temperature;

double density_full(double pR, double pTheta, double R, double Theta, double phi, double pPhi)
{
    double H = pR*pR/2.0/mu + pTheta*pTheta/2.0/mu/R/R + pPhi*pPhi/2.0/mu/R/R/std::sin(Theta)/std::sin(Theta) + \
               ar_he_pot(R);

    if ( H > 0 )
        return std::exp(-H/3.166812E-6/Temperature);

    else
        return 0.0;
}

/*
void syst_2D (REAL t ,REAL *y, REAL *f)
{
   (void)(t); // avoid unused parameter warning

  double * out = new double[4];
  rhs_2D(out, y[0],y[1],y[2],y[3] );// R, pR, psi, ppsi

  f[0] = out[0];//dR/dt
  f[1] = out[1];
  f[2] = out[2];//d(psi)/dt
  f[3] = out[3];//d(ppsi)/dt

  delete [] out;
}
*/

void syst_3D( REAL t, REAL * y, REAL * f )
{
    (void)(t);

    double * out = new double [6];
    rhs_3D(out, y[0], y[1], y[2], y[3], y[4], y[5]);

    f[0] = out[0]; // dR / dt
    f[1] = out[1]; // dpR / dt
    f[2] = out[2]; // d theta / dt
    f[3] = out[3]; // d pTheta / dt
    f[4] = out[4]; // d phi / dt
    f[5] = out[5]; // d pPhi / dt

    delete [] out;
}

void create_frequencies_vector( std::vector<double> & freqs, const double sampling_time )
{
    const double FREQ_MAX = 1000.0;
    double FREQ_STEP = 1.0 / (sampling_time * constants::ATU) / constants::LIGHTSPEED_CM / MaxTrajectoryLength; // cm^-1
    int FREQ_SIZE = (int) FREQ_MAX / FREQ_STEP + 1;

    freqs.resize( FREQ_SIZE );
    for( int k = 0; k <  FREQ_SIZE; ++k )
        freqs[k] = k * FREQ_STEP;
}

double spectral_function( const double pR, const double theta, const double pTheta,
                          const double phi, const double pPhi, const int freq_bin )
{
    REAL     epsabs = 1.0E-13;       /* absolute error bound                      */
    REAL     epsrel = 1.0E-13;       /* relative error bound                      */
    REAL     t0 = 0.0;           /* left edge of integration interval         */
    REAL     *y0;          /* [0..n-1]-vector: initial value, approxim. */

    REAL     h = 0.1;            /* initial, final step size                  */
    REAL     xend = sampling_time;         /* right edge of integration interval        */
    long     fmax = 1e8;         /* maximal number of calls of right side     */
                             /* in gear4()                                */
    long     aufrufe;      /* actual number of function calls           */

    int      N = 6;             //number of DEs in system
    int      fehler;       /* error code from umleiten(), gear4()       */

    void     *vmblock = vminit();     /* List of dynamically allocated vectors     */

    y0  = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);

    if (! vmcomplete(vmblock))  {       /* out of memory? */
        printf("mgear: out of memory.\n");
        throw std::runtime_error("out of memory!");
    }

    std::vector<double> dipx, dipy, dipz;
    double dipx_, dipy_, dipz_;

    y0[0] = R0_max * 0.99;
    y0[1] = pR;
    y0[2] = theta;
    y0[3] = pTheta;
    y0[4] = phi;
    y0[5] = pPhi;

    std::cout << "initial values: " << y0[0] << " " << y0[1] << " " << y0[2] << " " << y0[3]
              << " " << y0[4] << " " << y0[5] << std::endl;

    for ( int counter = 0; y0[0] < R0_max; ++counter )
    {
        fehler = gear4(&t0, xend, N, syst_3D, y0, epsabs, epsrel, &h, fmax, &aufrufe);

        if (fehler != 0)
        {
            std::cout << "Gear4: error n " << 10 + fehler << std::endl;
            break;
        }

        if ( counter == MaxTrajectoryLength )
        {
            std::cout << std::endl << "Trajectory cut" << std::endl << std::endl;
            return 0.0;
        }

        dipx_ = ar_he_dip_buryak_fit(y0[0]) * std::sin(y0[2]) * std::cos(y0[4]);
        dipy_ = ar_he_dip_buryak_fit(y0[0]) * std::sin(y0[2]) * std::sin(y0[4]);
        dipz_ = ar_he_dip_buryak_fit(y0[0]) * std::cos(y0[2]);

        dipx.push_back( dipx_ );
        dipy.push_back( dipy_ );
        dipz.push_back( dipz_ );

        xend = sampling_time * (counter + 2);
        aufrufe = 0;
    }

    std::cout << "trajectory length: " << dipx.size() << std::endl;

    double * dipx_in = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );
    double * dipy_in = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );
    double * dipz_in = (double*) fftw_malloc( sizeof(double) * MaxTrajectoryLength );

    std::fill( dipx_in, dipx_in + MaxTrajectoryLength, 0.0 );
    std::fill( dipy_in, dipy_in + MaxTrajectoryLength, 0.0 );
    std::fill( dipz_in, dipz_in + MaxTrajectoryLength, 0.0 );

    for ( size_t k = 0; k < dipx.size(); ++k )
    {
        dipx_in[k] = dipx[k];
        dipy_in[k] = dipy[k];
        dipz_in[k] = dipz[k];
    }

    fftw_complex * dipx_out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * FFTMaxLength );
    fftw_complex * dipy_out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * FFTMaxLength );
    fftw_complex * dipz_out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * FFTMaxLength );

    for ( int k = 0; k < FFTMaxLength; ++k )
    {
        dipx_out[k][0] = 0.0;
        dipx_out[k][1] = 0.0;

        dipy_out[k][0] = 0.0;
        dipy_out[k][1] = 0.0;

        dipz_out[k][0] = 0.0;
        dipz_out[k][1] = 0.0;
    }

    fftw_plan px = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, dipx_in, dipx_out, FFTW_ESTIMATE );
    fftw_plan py = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, dipy_in, dipy_out, FFTW_ESTIMATE );
    fftw_plan pz = fftw_plan_dft_r2c_1d( MaxTrajectoryLength, dipz_in, dipz_out, FFTW_ESTIMATE );

    fftw_execute( px );
    fftw_execute( py );
    fftw_execute( pz );

    //double specfunc_coeff = 1.0/(4.0*M_PI)/constants::EPSILON0 * std::pow(sampling_time * constants::ATU, 2)/2.0/M_PI * \
    //        std::pow(constants::ADIPMOMU, 2);

    double dipfft = dipx_out[freq_bin][0] * dipx_out[freq_bin][0] + dipx_out[freq_bin][1] * dipx_out[freq_bin][1] + \
                    dipy_out[freq_bin][0] * dipy_out[freq_bin][0] + dipy_out[freq_bin][1] * dipy_out[freq_bin][1] + \
                    dipz_out[freq_bin][0] * dipz_out[freq_bin][0] + dipz_out[freq_bin][1] * dipz_out[freq_bin][1];

    // j -> erg -- 1E7; m3 -> cm6 -- 1E12
    //double specfunc_value = 1.0E19 * specfunc_coeff * dipfft;
    //std::cout << "specfunc_value: " << specfunc_value << std::endl;
    double specfunc_value = dipfft;

    fftw_destroy_plan( px );
    fftw_destroy_plan( py );
    fftw_destroy_plan( pz );


    fftw_free( dipx_in );
    fftw_free( dipy_in );
    fftw_free( dipz_in );
    fftw_free( dipx_out );
    fftw_free( dipy_out );
    fftw_free( dipz_out );

    vmfree( vmblock );

    return specfunc_value;
}

double integrand_full( hep::mc_point<double> const& x, const int freq_bin)
{
    double pR = std::tan(M_PI * (x.point()[0] - 0.5));
    double Theta = x.point()[1] * M_PI;
    double pTheta =  std::tan(M_PI * (x.point()[2] - 0.5));
    double phi = x.point()[3] * 2.0 * M_PI;
    double pPhi = std::tan(M_PI * (x.point()[4] - 0.5));

    double pRJ = M_PI * (1.0 + pR*pR);
    double THJ = M_PI;
    double pTJ = M_PI * (1 + pTheta*pTheta);
    double PHJ = 2.0 * M_PI;
    double PPHJ = M_PI * (1.0 + pPhi*pPhi);

    double result = std::abs(pR) / mu * spectral_function(pR, Theta, pTheta, phi, pPhi, freq_bin ) * \
            density_full(pR, pTheta, R0_max, Theta, phi, pPhi);

    return result * THJ * pRJ * pTJ * PHJ * PPHJ;
}

double denumerator_integrand( hep::mc_point<double> const& x )
{
    double R = (R0_max - R0_min) * x.point()[0];
    double pR = std::tan( M_PI * (x.point()[1] - 0.5 ));
    double Theta = M_PI * x.point()[2];
    double pTheta = std::tan( M_PI * (x.point()[3] - 0.5) );
    double phi = 2.0 * M_PI * x.point()[4];
    double pPhi = std::tan( M_PI * (x.point()[5] - 0.5) );

    double RJ = (R0_max - R0_min);
    double PRJ = M_PI * (1.0 + pR * pR);
    double THJ = M_PI;
    double PTJ = M_PI * (1.0 + pTheta * pTheta);
    double PHJ = 2.0 * M_PI;
    double PPHJ = M_PI * (1.0 + pPhi * pPhi);

    double result = density_full(pR, pTheta, R, Theta, phi, pPhi);

    return result * RJ * PRJ * THJ * PTJ * PHJ * PPHJ;
}

int main( int argc, char* argv[] )
{
    // initialize MPI
    MPI_Init( &argc, &argv );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

    std::ofstream file;
    file.open( "spectral_function.txt" );

    hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

    auto denominator_results = hep::mpi_vegas(
         MPI_COMM_WORLD,
         hep::make_integrand<double>(denumerator_integrand, 6),
         std::vector<std::size_t>(10, 1e4)
    );

    auto denominator_result = hep::cumulative_result0( denominator_results.begin() + 2, denominator_results.end() );

    std::cout << "denominator result: " << denominator_result.value() << std::endl;

    std::vector<double> freqs;
    create_frequencies_vector( freqs, sampling_time );

    double MINFREQ = 0.0;
    double MAXFREQ = 0.0;
    double freqstep = 5.0;

    for ( double FREQ = MINFREQ; FREQ <= MAXFREQ; FREQ += freqstep )
    {
        std::vector<double>::iterator it = std::lower_bound( freqs.begin(), freqs.end(), FREQ );
        int freq_bin = (int) (it - freqs.begin());

        if ( rank == 0 )
        {
            std::cout << "-------------- spectral function calculation; freq = "<< FREQ << " ----------" << std::endl;
            std::cout << "FREQ: " << FREQ << "; freq_bin: " << freq_bin << "; freq_bin value: " << freqs[freq_bin] << std::endl;
        }

        std::chrono::high_resolution_clock::time_point cycleStartTime = std::chrono::high_resolution_clock::now();

        auto numerator_integrand = std::bind( integrand_full, std::placeholders::_1, freq_bin );

        auto numerator_results = hep::mpi_vegas(
            MPI_COMM_WORLD,
            hep::make_integrand<double>(numerator_integrand, 5),
            std::vector<std::size_t>(5, 4e2)
        );

        auto numerator_result = hep::cumulative_result0( numerator_results.begin() + 2, numerator_results.end() );

        //double VOLUME = 4.0/3.0 * M_PI * std::pow(R0_max * constants::ALU, 3.0);
        double VOLUME = 4.0/3.0 * M_PI * std::pow(R0_max, 3.0);
        double result = VOLUME * numerator_result.value() / denominator_result.value();

        // в лоб переводим размерность итоговой величины
        result = result * constants::HTOJ * std::pow(constants::ALU, 6) * constants::ATU * 1.0E19;

        if ( rank == 0 )
        {
            file <<  FREQ << " " <<  std::setprecision(10) << result << std::endl;
            std::cout << "Final result: " << result<< std::endl;
            std::chrono::high_resolution_clock::time_point cycleEndTime = std::chrono::high_resolution_clock::now();

            std::cout << "Cycle time: " <<
                std::chrono::duration_cast<std::chrono::milliseconds>(cycleEndTime - cycleStartTime).count() / 1000.0
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
