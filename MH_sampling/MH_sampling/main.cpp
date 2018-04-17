#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>

#include "ar_he_dip_buryak_fit.hpp"
#include "ar_he_pes.hpp"

const double Racc_max = 60.0;
const double Racc_min = 2.0;
const int dim = 6;

const double Rint_max = 40.0;
const double Rint_min = 4.0;

const double a0 = 0.529E-10; //in m

const double VOLUME = 4.0/3.0*M_PI*pow(40.0*a0,3.0)  ;

const double mu = 6632;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 generator{ seed };

const double r_step = 0.6;
const double angles_step = 0.6;
const double impulse_step = 10.0;


//x[0] == pR, x[1] == pTheta, x[2] = R, x[3] == pPhi, x[4] == Theta, x[5] == Phi
void sample( std::vector<double> const & prev, std::vector<double> & next )
{
    // Different gaussians in each direction
    /*
    std::normal_distribution<double> distribution_r(0.0, r_step);
    std::normal_distribution<double> distribution_angles(0.0, angles_step);
    std::normal_distribution<double> distribution_impulses(0.0, impulse_step);
    next[0] = prev[0] + distribution_impulses( generator );
    next[1] = prev[1] + distribution_impulses( generator );
    next[2] = prev[2] + distribution_r( generator );
    next[3] = prev[3] + distribution_impulses( generator );
    next[4] = prev[4] + distribution_angles( generator );
    next[5] = prev[5] + distribution_angles( generator );
    */
    // ----------------------------------------------------
    // ATMCMC -- https://arxiv.org/pdf/1408.6667.pdf ФИГНЯ?
    /*
    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    // proposal distribution
    std::normal_distribution<double> q(0.0, 2.0);
    double epsilon = q(generator);

    for ( int i = 0; i < dim; ++i )
        next[i] = prev[i] + epsilon * unidistr(generator);
    */
    // ----------------------------------------------------

    // Same gaussian for each direction
    // std::normal_distribution<double> distribution_gen(0.0, 1.5);
    std::normal_distribution<double> distribution_gen(0.0, 0.5);
    for (int i = 0; i < dim; ++i)
        next[i] = prev[i] + distribution_gen( generator );
}

double integrand(const double x)
{
    return  ar_he_dip_buryak_fit(x) * ar_he_dip_buryak_fit(x);
}

// W = \rho * |mu|^2
//x[0] == pR, x[1] == pTheta, x[2] = R, x[3] == pPhi, x[4] == Theta, x[5] == Phi
double density( std::vector<double> const & x )
{
    double H = x[0] * x[0] / 2.0 / mu + \
            x[1] *x[1] / mu / x[2] / x[2] / 2.0 + \
            x[3] * x[3] / mu / x[2] / x[2] / std::sin(x[4]) / std::sin(x[4])
     + ar_he_pot(x[2]);

    if ((H > 0)  && (x[2] > Racc_min) && (x[2] < Racc_max))
    {
        return ar_he_dip_buryak_fit(x[2]) * ar_he_dip_buryak_fit(x[2]) * exp(-H/3.166812E-6/295.0);
    }
    else
        return 0.0;
}

double map_ang(double ang, double lim)
{
    return fmod(lim + fmod(ang,lim), lim);
}

double MHA(const int burnin, const int chain_length, std::vector<double> & curr, std::string const & filename )
{
    std::vector<double> xcand(dim);

    double alpha, u;
    int accepted = 0;
    int integral_counter = 0;

    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    // ----------------------------------
    // Burnin
    for (int i = 0; i < burnin; ++i)
    {
        sample(curr, xcand);
        alpha = density(xcand) / density(curr);
        u =  unidistr( generator );

        if ( u < alpha )
            curr = xcand;
    }
    // ----------------------------------

    double integral = 0.0;
    double mean = 0.0, d;

    std::ofstream out(filename);

    int iterations = 0;
    for ( ; integral_counter < chain_length ; ++iterations )
    {
        sample(curr, xcand);
        alpha = density(xcand) / density(curr);
        u =  unidistr( generator );
        if ( u < std::min(alpha,1.0) )
        {
            curr = xcand;
            ++accepted;
        }

        if (iterations % 20 == 0)
        {
            if ( (curr[2] > Rint_min) && (curr[2] < Rint_max) )
            {
                 out << curr[2] << " " << curr[0] << " " << map_ang(curr[4], M_PI) << " " << curr[1] << " "
                     << map_ang(curr[5], 2 * M_PI) << " " << curr[3] << std::endl;

                 integral += integrand(curr[2]);
                 d = integrand(curr[2]) - mean;
                 mean += d / ( (double) integral_counter + 1.0 );
                 ++integral_counter;

                 if ( integral_counter % 1000 == 0 )
                    std::cout << "integral_counter: " << integral_counter << std::endl;
            }
        }

        //if ( iterations % 10000 == 0 )
        //    std::cout << ">> iterations = " << iterations << "; integral_counter = " << integral_counter << "..." << std::endl;
    }

    double integral_value = integral/integral_counter; // * VOLUME;

    std::cout << "Percentage of accepted: " << accepted / (double) iterations * 100.0 << std::endl;
    // std::cout << "Dipole mean over W: " << integral_value << std::endl;

    return integral_value;
}

int main()
{
    std::vector<double> x0{0.0, 0.0, 7.0, 0.0, 0.0, 0.0};

    const std::string filename = "../initial_points_zimm_50000.txt";
    const int burnin = 1e5;
    const int chain_len = 5e4;

    auto start = std::clock();

    const int ncycles = 1;
    std::vector<double> vals(ncycles);
    for ( int i = 0; i < ncycles; ++i )
        vals[i] = MHA(burnin, chain_len, x0, filename);

    double mean =std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    std::cout << "Mean value: " << mean << std::endl;

    double sq_sum = std::inner_product(vals.begin(), vals.end(), vals.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / vals.size() - mean * mean);
    std::cout << "Standard deviation: " << stdev << std::endl;

    std::cout << "Time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    return 0;
}

