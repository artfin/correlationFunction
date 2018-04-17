#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <utility> // std::pair, std::make_pair
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include "constants.hpp"

const double Temperature = 295.0;

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

double integrate_simpson(double omega, std::vector<std::pair<double, double>> const & input )
{
    double integral = 0.0;
    double step = (input.end()[-1].first - input[0].first) / input.size();

    for (size_t i = 1; i < input.size() - 1; i += 2)
    {
        integral += integrand(input[i - 1].first, input[i - 1].second, omega) +
                    4 * integrand(input[i].first, input[i].second, omega) +
                    integrand(input[i + 1].first, input[i + 1].second, omega);

    }
    std::cout << "integral: " << integral << std::endl;

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

// симметризация корреляционной функции и частичное преобразование размерности
void symmetrize( std::vector<std::pair<double, double>> const & input, std::vector<std::pair<double, double>> & res )
{
    res.resize( 2 * input.size() - 1 );
    std::cout << "input.size(): " << input.size() << std::endl;
    std::cout << "res.size(): " << res.size() << std::endl;

    int j = 0;
    for ( int i = (int) input.size() - 1; i >= 0; --i, ++j )
        res[j] = std::make_pair( -input[i].first, constants::ADIPMOMU * constants::ADIPMOMU / \
                                 (4.0 * M_PI * constants::EPSILON0) * input[i].second );

    for ( int i = 1; i < (int) input.size(); ++i, ++j )
        res[j] = std::make_pair( input[i].first, constants::ADIPMOMU * constants::ADIPMOMU / \
                                 (4.0 * M_PI * constants::EPSILON0) * input[i].second );
}

int main( )
{
    std::vector<std::pair<double, double>> input;
    std::ifstream inFile("../data/input.txt");
    reader(inFile, input);

    // симметризованная корреляционная функция
    // размерность: Дж * м^6
    std::vector<std::pair<double, double>> symm_input;
    symmetrize( input, symm_input );

    std::ofstream symFile ("../data/symm.txt");
    for ( int i = 0; i < (int) symm_input.size(); ++i )
        symFile << symm_input[i].first << " " << symm_input[i].second << std::endl;
    symFile.close();

    int Nomega = 1000;
    double omega0  = 0.0;
    double omegastep = 1500.0/Nomega;

    std::ofstream outFile("../data/output.txt");
    std::ofstream specFile("../data/spectrum.txt");
    for (int i = 0; i < Nomega; ++i)
    {
        double specfunc = integrate_simpson(omega0, input);
        outFile << omega0 << " " << specfunc << std::endl;

        double spectrum = std::pow(2*M_PI,3.0)*constants::LOSHMIDT_CONSTANT*constants::LOSHMIDT_CONSTANT/3.0/ \
                constants::PLANCKCONST*omega0*(1-exp(-2*0.7193875*omega0/295.0))*specfunc;

        //double spectrum = std::pow(2.0 * M_PI, 3.0) * std::pow(constants::LOSHMIDT_CONSTANT, 2) / 3.0 / constants::PLANCKCONST_REDUCED \
        //        * std::pow(10.0, -2.0) * omega0 * (1 - std::exp(- 2 * 0.7193875 * omega0 / Temperature)) * specfunc;

        // 1/c * (c*omega0) = 10**(-2) * omega0
        // omega0 -- в cm^-1, а для соблюдения размерности формулы должна быть в Hz => переводной коэффициент "c" сокращается

        specFile << omega0 << " " << spectrum << std::endl;

        omega0 += omegastep;
    }

    outFile.close();
    specFile.close();

    return 0;
}

