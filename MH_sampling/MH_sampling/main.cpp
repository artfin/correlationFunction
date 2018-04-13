#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <random>
#include <algorithm>
#include <vector>
#include <chrono>

#include "ar_he_dip_buryak_fit.hpp"
#include "ar_he_pes.hpp"

using std:: uniform_real_distribution;

double step = 1.5;

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

void sample( std::vector<double> const & prev, std::vector<double> & next )
{
    std::normal_distribution<double> distribution_gen( 0.0, step );
    for (int i = 0; i < dim; ++i)
         next[i] = prev[i] + distribution_gen( generator );
}

double integrand(const double x)
{
    return  ar_he_dip_buryak_fit(x) * ar_he_dip_buryak_fit(x);
}

//x[0] == pR, x[1] == pTheta, x[2] = R, x[3] == pPhi, x[4] == Theta, x[5] == Phi
double density( std::vector<double> const & x )
{
    double H = x[0] * x[0] / 2.0 / mu + \
            x[1] *x[1] / mu / x[2] / x[2] / 2.0 + \
            x[3] * x[3] / mu / x[2] / x[2] / std::sin(x[4]) / std::sin(x[4])
     + ar_he_pot(x[2]);

    if ((H > 0)  && (x[2] > Racc_min) && (x[2] < Racc_max))
    {
        return exp(-H/3.166812E-6/295.0);
    }
    else
        return 0.0;
}

double map_ang(double ang, double lim)
{
    return fmod(lim + fmod(ang,lim), lim);
}

double MHA(int burnin,int iterations, FILE * flname, std::vector<double> & curr )
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

    for ( int counter = 0; counter < iterations; ++counter )
    {
        sample(curr, xcand);
        alpha = density(xcand) / density(curr);
        u =  unidistr( generator );
        if ( u < std::min(alpha,1.0) )
        {
            curr = xcand;
            ++accepted;
        }

        if (counter % 10 == 0)
        {
            //fprintf(flname, "%f\n", map_ang(xi_1[5], 2*M_PI));
            //fprintf(flname, "%f\n", (xi_1[2]));
            if ( (curr[2] > Rint_min) && (curr[2] < Rint_max) )
            {
                 integral += integrand(curr[2]);
                 d = integrand(curr[2]) - mean;
                 mean += d / ( (double) integral_counter + 1.0 );
                 ++integral_counter;
            }
        }
    }

    double integral_value = integral/integral_counter * VOLUME;
    printf("The percentage of accepted: %f\n", accepted / (double)iterations * 100.0 );
    printf("Integration: %e %e\n", integral_value, mean * VOLUME);

    return integral_value;
}





int main()
{
    FILE * f1 = fopen("weight_density.txt", "w");

    std::vector<double> x0{0.0, 0.0, 7.0, 0.0, 0.0, 0.0};

    double mean = 0.0;
    int cycles = 5;

    for (int i = 0; i < cycles; ++i)
        mean  += MHA(30000,5000000, f1, x0);

    printf("mean: %e\n", mean/(double)cycles);
    fclose(f1);

    return 0;
}

