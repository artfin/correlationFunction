#include <iostream>
#include <random>
#include <gsl/gsl_histogram.h>
#include <fstream>
#include <functional>
#include <chrono>
#include <vector>
#include "ar_he_pot.hpp"
#include "constants.hpp"

unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 mcmc_generator{ seed };

const double step = 1.5;
const double DIM = 6;

const double Racc_max = 60.0;
const double Racc_min = 2.0;

const double Rint_max = 40.0;
const double Rint_min = 4.0;

const double MU = 6632.039;

void sample(std::vector<double> const & prev, std::vector<double> & next)
{
    std::normal_distribution<double> distribution( 0.0, step );
    for ( size_t i = 0; i < next.size(); ++i )
         next[i] = prev[i] + distribution( mcmc_generator );
}

// порядок перменных: R, pR, theta, pT, phi, pPhi
double density_( std::vector<double> & x, const double Temperature )
{
    double R = x[0];
    double pR = x[1];
    double theta = x[2];
    double pTheta = x[3];
    // double phi = x[4];
    double pPhi = x[5];

    double H = std::pow(pR, 2) / (2.0 * MU) + \
               std::pow(pTheta, 2) / (2.0 * MU * R * R) + \
               std::pow(pPhi, 2) / (2.0 * MU * R * R * sin(theta) * sin(theta)) + \
               ar_he_pot(R);

    if ( (H > 0) && (R > Racc_min) && (R < Racc_max) )
        return exp(- H * constants::HTOJ / (constants::BOLTZCONST * Temperature) );
    else
        return 0.0;
}

std::vector<double> MHA_burnin( int burnin, std::vector<double> x, std::function<double(std::vector<double>&)> density )
{
    std::vector<double> xcand(DIM);
    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    for (int i = 0; i < burnin; ++i)
    {
        sample(x, xcand);
        double alpha = density(xcand) / density(x);

        if ( unidistr(mcmc_generator) < alpha)
        {
            x = xcand;
        }
    }

    return x;
}

// getOutOf -- то количество точек, раз в которое мы выбираем точку,
// чтобы на ней посчитать интегранд
std::vector<double> MHA_generate_point(std::vector<double> x, int getOutOf,
                                       std::function<double(std::vector<double>&)> density)
{
    std::vector<double> xcand(DIM);
    int counter = 1;

    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    while( true)
    {
        sample(x, xcand);
        double alpha = density(xcand) / density(x);

        if ( unidistr(mcmc_generator) < alpha )
        {
            x = xcand;
        }

        if ( (counter > getOutOf) && (xcand[0] > Rint_min) && (xcand[0] < Rint_max) )
        {
            return xcand;
        }

        ++counter;
    }
}

int main()
{
    std::cout << "---- Initializing Metropolis-Hastings sampler -------------" << std::endl;

    std::function<double(std::vector<double>&)> density = bind( density_, std::placeholders::_1, 295.0 );

    const int burnin = 3e4;
    std::vector<double> initial_point = {6.0, 0.0, 1.2, 0.0, 0.0, 0.0 };
    std::vector<double> initial_after_burnin = MHA_burnin( burnin, initial_point, density );

    const int getOutOf = 10;
    std::vector<double> p = MHA_generate_point( initial_after_burnin, getOutOf, density );

    const int chainLength = 1e4;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for ( int k = 0; k < chainLength; k++ )
    {
       p = MHA_generate_point( p, getOutOf, density );
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double) 1000.0
              << " s" << std::endl;
    return 0;
}
