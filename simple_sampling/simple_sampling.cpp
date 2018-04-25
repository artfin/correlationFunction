#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <chrono>
#include <cmath>

#include "ar_he_pes.hpp"
#include "ar_he_dip_buryak_fit.hpp"
#include "constants.hpp"

const double Temperature = 295.0;
const int DIM = 6;
const double MU = 6632.039;

class Sampler
{
public:
    void init( const double Rmin, const double Rmax )
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator = std::mt19937(seed);

        //double Rmin_transformed = 1.0 / Rmin / Rmin;
        //double Rmax_transformed = 1.0 / Rmax / Rmax;
        RDistr = std::uniform_real_distribution<double>( Rmin, Rmax );

        uniDistr = std::uniform_real_distribution<double>( -1.0, 1.0 );

        thetaDistr = std::uniform_real_distribution<double>( 0.0, M_PI );
        phiDistr = std::uniform_real_distribution<double>( 0.0, 2 * M_PI );
    }

    // R pR theta pT phi pPhi
    void sample( std::vector<double> & x )
    {
        x[0] = RDistr( generator );
        x[1] = uniDistr( generator );
        x[2] = thetaDistr( generator );
        x[3] = uniDistr( generator );
        x[4] = phiDistr( generator );
        x[5] = uniDistr( generator );
    }

private:
	std::uniform_real_distribution<double> RDistr;
	std::uniform_real_distribution<double> thetaDistr;
    std::uniform_real_distribution<double> phiDistr;

    std::uniform_real_distribution<double> uniDistr;

    std::mt19937 generator;
};

// R pR theta pTheta phi pPhi
double integrand( std::vector<double> const & x )
{
    double R = x[0];
    double pR = std::tan( M_PI * (x[1] - 0.5));
    double pTheta = std::tan( M_PI * (x[3] - 0.5) );
    double pPhi = std::tan( M_PI * (x[5] - 0.5) );

    double h = pR * pR / (2.0 * MU) + pTheta * pTheta / (2.0 * MU * R * R) + \
            pPhi * pPhi / (2.0 * MU * R * R * std::sin(x[2]) * std::sin(x[2])) + ar_he_pot(R);

    if ( h > 0 )
    {
        double PRJ = M_PI * (1.0 + pR * pR);
        double PTJ = M_PI * (1.0 + pTheta * pTheta);
        double PPJ = M_PI * (1.0 + pPhi * pPhi);
        double RJ = 1.0; //2.0 * std::pow(R, -3.0);

        double jacs = PRJ * PTJ * PPJ * RJ;
        //std::cout << "jacs: " << jacs << std::endl;
        //std::cout << "h: " << h << std::endl;

        double dip = ar_he_dip_buryak_fit(R);

        double res = jacs * dip * dip * std::exp(-h * constants::HTOJ / (constants::BOLTZCONST * Temperature));
        //std::cout << "res: " << res << std::endl;

        return res;
    }
    else
        return 0.0;
}


int main()
{
	std::cout << "-------------------------------------------" << std::endl;
	std::cout << "Program for generating uniform initial conditions for diatom system" << std::endl;

    const double Rmin = 6.0;
    const double Rmax = 8.0;

	std::cout << "Rmin: " << Rmin << "; Rmax: " << Rmax << std::endl;

    std::cout << "Tranforming pR, pTheta, pPhi to [-1, 1]" << std::endl;

    Sampler sampler;
    sampler.init( Rmin, Rmax );

    std::cout << "Testing sampling procedure...." << std::endl;
    std::cout << "Integrating mu^2 exp(-H/kT)" << std::endl;

    const int ncycles = 5;
    const int neval = 1e6;
    std::cout << "Running ncycles = " << ncycles << " using neval = " << neval << " function evaluations." << std::endl << std::endl;

    std::vector<double> results;
    std::vector<double> x(DIM);

    for ( int ncycle = 0; ncycle < ncycles; ++ncycle )
    {
        std::cout << ">> Starting iteration " << ncycle << std::endl;
        double sum = 0.0;

        for ( int iter = 0; iter < neval; ++iter )
        {
            sampler.sample(x);
            sum += integrand( x );
        }

        sum /= neval;
        results.push_back( sum );

        std::cout << ">> Integral value = " << sum << std::endl << std::endl;
    }

	return 0;
}
