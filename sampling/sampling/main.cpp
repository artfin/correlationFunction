#include <iostream>
#include <random>
#include <gsl/gsl_histogram.h>
#include <fstream>
#include <functional>
#include <chrono>
#include <vector>
#include <cassert>
#include "ar_he_pot.hpp"
#include "constants.hpp"

unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 mcmc_generator{ seed };

const double step = 2.5;
const double DIM = 6;
const int NBINS = 100;

const double Racc_max = 41.0;
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

    // !!!!!!!!!!!!!!!!!!!!!!!!!
    if ( (H < 0) && (R > Racc_min) && (R < Racc_max) )
        return exp(- H * constants::HTOJ / (constants::BOLTZCONST * Temperature) );
    else
        return 0.0;
}

std::vector<double> MHA_burnin( int burnin, std::vector<double> x, std::function<double(std::vector<double>&)> density )
{
    std::vector<double> xcand(DIM);
    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    int accepted = 0;
    for (int i = 0; i < burnin; ++i)
    {
        sample(x, xcand);
        double alpha = density(xcand) / density(x);

        if ( unidistr(mcmc_generator) < alpha)
        {
            x = xcand;
            ++accepted;
        }
    }

    std::cout << "(burnin) Acceptance rate: " << (double) accepted / burnin * 100 << "%" << std::endl;

    return x;
}

double wrapMax( double x, double xmax )
{
    return std::fmod( xmax + fmod(x, xmax), xmax );
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
            // r pr theta[2] ptheta varphi[4] pvarphi
            xcand[2] = wrapMax(xcand[2], M_PI);
            xcand[4] = wrapMax(xcand[4], 2 * M_PI);
            return xcand;
        }

        ++counter;
    }
}

void allocateHistograms( std::vector<gsl_histogram*> & histograms,
                         std::vector<std::pair<double, double>> & ranges, const int nbins )
{
    for ( int k = 0; k < DIM; k++ )
    {
        gsl_histogram * h = gsl_histogram_alloc(nbins);
        gsl_histogram_set_ranges_uniform(h, ranges[k].first, ranges[k].second);
        histograms.push_back( h );
    }
}

void addPoint( std::vector<gsl_histogram*> & histograms, std::vector<double> & p )
{
    assert( p.size() == histograms.size() );
    for ( size_t k = 0; k < histograms.size(); k++ )
        gsl_histogram_increment(histograms[k], p[k]);
}

void normalizeHistogram( gsl_histogram * h )
{
    double sum = 0.0;
    double lower_bound, upper_bound;
    for ( size_t k = 0; k < gsl_histogram_bins(h); k++ )
    {
        gsl_histogram_get_range(h, k, &lower_bound, &upper_bound);
        double step = upper_bound - lower_bound;
        double content = gsl_histogram_get(h, k);
        sum += content * step;
    }

    gsl_histogram_scale(h, 1.0 / sum);
}

void normalizeHistograms( std::vector<gsl_histogram*> histograms )
{
    for ( size_t k = 0; k < histograms.size(); k++ )
        normalizeHistogram(histograms[k]);
}

void saveHistogram( gsl_histogram * h, std::string filename )
{
    std::ofstream outFile(filename);
    int bins = gsl_histogram_bins(h);
    double lower_bound, upper_bound;
    for ( int k = 0; k < bins; k++ )
    {
        gsl_histogram_get_range(h, k, &lower_bound, &upper_bound);
        double bin_center = 0.5 * (lower_bound + upper_bound);
        double content = gsl_histogram_get(h, k);
        outFile << bin_center << " " << content << std::endl;
    }

    outFile.close();
}

void saveHistograms( std::vector<gsl_histogram*> histograms )
{
    std::cout << ">> Saving histograms..." << std::endl;
    std::vector<std::string> names = {"r.txt", "pr.txt", "theta.txt", "ptheta.txt", "varphi.txt", "pvarphi.txt"};
    assert( histograms.size() == names.size() );

    std::string folder = "../hist_bound/";

    for ( size_t k = 0; k < names.size(); k++ )
        saveHistogram(histograms[k], folder + names[k]);
}

void freeHistograms( std::vector<gsl_histogram*> histograms )
{
    std::cout << ">> Freeing histograms..." << std::endl;
    for ( size_t k = 0; k < histograms.size(); ++k )
        gsl_histogram_free(histograms[k]);
}

void show_progress( double progress )
{
    std::cout << "[";
    int barWidth = 70;
    int pos = barWidth * progress;
    for ( int i = 0; i < barWidth; ++i )
    {
        if ( i < pos ) std::cout << "=";
        else if ( i == pos ) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

int main()
{
    std::cout << "---- Initializing Metropolis-Hastings sampler -------------" << std::endl;

    std::function<double(std::vector<double>&)> density = bind( density_, std::placeholders::_1, 295.0 );

    const int burnin = 3e4;
    std::vector<double> initial_point = {6.0, 0.0, 1.2, 0.0, 0.0, 0.0 };
    std::vector<double> initial_after_burnin = MHA_burnin( burnin, initial_point, density );
    std::cout << ">> Burnin is done." << std::endl;

    const int getOutOf = 1e3;
    std::vector<double> p = MHA_generate_point( initial_after_burnin, getOutOf, density );

    // аллокируем гистограммы
    std::vector<gsl_histogram*> histograms;
    std::vector<std::pair<double, double>> ranges;
    ranges.emplace_back(4.0, 40.0); // R
    ranges.emplace_back(-12.0, 12.0); // pR
    ranges.emplace_back(0.0, M_PI); // theta
    ranges.emplace_back(-250.0, 250.0); // pTheta
    ranges.emplace_back(0.0, 2 * M_PI); // varphi
    ranges.emplace_back(-250.0, 250.0); // pVarphi
    allocateHistograms( histograms, ranges, NBINS );

    const int chainLength = 10;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    std::ofstream out("../samples_bound.txt");

    const int blockSize = 1e3;
    std::cout << "[";
    for ( int k = 0; k < chainLength; k++ )
    {
        p = MHA_generate_point( p, getOutOf, density );
        addPoint(histograms, p);

        for ( int k = 0; k < DIM; ++k )
            out << p[k] << " ";
        out << std::endl;

        if ( (k % blockSize == 0) == 0 )
            show_progress((double) k / chainLength);
    }

    out.close();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double) 1000.0
              << " s" << std::endl;

    normalizeHistograms(histograms);
    saveHistograms(histograms);

    std::cout << "---- Quitting Metropolis-Hastings sampler --------" << std::endl;

    return 0;
}
