#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include "constants.hpp"

void reader( std::ifstream & inFile, std::vector<std::pair<double, double>> & contents )
{
    if ( !inFile )
        std::invalid_argument("Can't open the file");

    const int MAXLINE = 100;
    char buffer[MAXLINE];
    std::stringstream ss;
    double time, corr;

    while( inFile.getline(buffer, MAXLINE ))
    {
        ss << buffer;
        ss >> time >> corr;
        contents.emplace_back(time, corr);
        ss.clear();
        ss.str("");
    }
}

double integrand(double time, double func, double omega)
{
    return func * std::cos(time * constants::LIGHTSPEED_CM * omega * 2 * M_PI);
}

double integrate(double omega, std::vector<std::pair<double, double>> const & input )
{
    double integral = 0.0;
    double step = (input.end()[-1].first - input[0].first) / input.size();

    for (size_t i = 1; i < input.size() - 1; i += 2)
    {
        integral += integrand(input[i - 1].first, input[i - 1].second, omega) +
                    4 * integrand(input[i].first, input[i].second, omega) +
                    integrand(input[i + 1].first, input[i + 1].second, omega);
    }

    integral *= step / 3.0;
    return integral;
}

double integrate_trapezoid( double omega, std::vector<std::pair<double, double>> & input )
{
    double integral = 0.0;
    double step = (input.end()[-1].first - input[0].first) / input.size();

    for ( size_t i = 1; i < input.size(); ++i )
    {
        integral += 0.5 * ( integrand(input[i].first, input[i].second, omega) + \
                            integrand(input[i - 1].first, input[i - 1].second, omega) );
    }

    integral *= step;
    return integral;
}

int main( )
{
    std::vector<std::pair<double, double>> input;
    std::ifstream inFile("../data/input.txt");
    reader(inFile, input);

    int Nomega = 1000;
    double omega0  = 0.0;
    double omegastep = 1500.0/Nomega;

    std::ofstream outFile("../data/output.txt");
    // outFile << std::fixed << std::setprecision(6);

    for (int i = 0; i < Nomega; ++i)
    {
        outFile << omega0 << " " << integrate_trapezoid(omega0, input) << std::endl;
        omega0 += omegastep;
    }

    return 0;
}

