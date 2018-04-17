#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <hep/mc.hpp>
#include "ar_he_pot.hpp"
#include "ar_he_dip_buryak_fit.hpp"

// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);
// avogadro number
const double AVOGADRO  = 6.022140857 * pow(10, 23);
// atomic length unit
const double ALU = 5.291772 * pow(10, -11);

const double MU = 6632.039;

// ------------------------------------------------------
// расстояние, вплоть до которого ведется интегрирование
const double Rint_min = 4.0; // a0
const double Rint_max = 25.0; // a0
// ------------------------------------------------------

double integrand_(hep::mc_point<double> const& x, double Temperature, bool numerator )
{
    double R_new = x.point()[0];
    double pR_new = x.point()[1];
    double theta_new = x.point()[2];
    double pTheta_new = x.point()[3];
    double pPhi_new = x.point()[5];

    double R = std::tan( M_PI / 2 * R_new );
    double Theta = M_PI * theta_new;
    double pR = std::tan(M_PI * (pR_new - 0.5));
    double pTheta = std::tan(M_PI * (pTheta_new - 0.5));
    double pPhi = std::tan(M_PI * (pPhi_new - 0.5));

    double potential_value = ar_he_pot(R);
    // std::cout << "R: " << R << "; pot: " << potential_value << std::endl;

    double h = std::pow(pR, 2) / (2 * MU) + std::pow(pTheta, 2) / (2.0 * MU * std::pow(R, 2)) + \
               std::pow(pPhi, 2.0) / (2.0 * MU * std::pow(R, 2.0) * std::pow(std::sin(Theta), 2.0)) + \
               potential_value;

    if ( R > Rint_max || R < Rint_min )
        return 0.0;


    if ( h > 0 )
    {
        double jacR = M_PI / 2.0 * (1.0 + std::pow(R, 2.0));
        double jacTheta = M_PI;
        double jacPhi = 2.0 * M_PI;
        double jacpR = M_PI * (1.0 + std::pow(pR, 2.0));
        double jacpTheta = M_PI * (1.0 + std::pow(pTheta, 2.0));
        double jacpPhi = M_PI * (1.0 + std::pow(pPhi, 2.0));

        // оба интеграла по области h > 0
        // в числителе стоит интеграл \mu^4 \exp(-H/kT)
        if ( numerator )
        {
            double dip = ar_he_dip_buryak_fit(R);
            double res = std::pow(dip, 4) * std::exp(-h * HTOJ / (BOLTZCONST * Temperature));
            return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
        }
        // в знаменателе стоит интеграл \mu^2 \exp(-H/kT)
        else
        {
            double dip = ar_he_dip_buryak_fit(R);
            double res = std::pow(dip, 2) * std::exp(-h * HTOJ / (BOLTZCONST * Temperature));
            return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
        }
    }

    else
        return 0.0;
}

int main()
{
    const double Temperature = 295.0;
    std::cout << "Current temperature = " << Temperature << std::endl;

    // set the verbose vegas callback function
    hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

    auto num_integrand = bind( integrand_, std::placeholders::_1, Temperature, true );
    auto denom_integrand = bind( integrand_, std::placeholders::_1, Temperature, false );

    auto results = hep::vegas(
       hep::make_integrand<double>(num_integrand, 6),
       std::vector<std::size_t>(10, 1e6)
     );

    auto result = hep::cumulative_result0(results.begin() + 1, results.end());
    double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());
    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;

    double numerator = result.value();

    results = hep::vegas(
         hep::make_integrand<double>(denom_integrand, 6),
         std::vector<std::size_t>(10, 1e6)
                );

    result = hep::cumulative_result0(results.begin() + 1, results.end());

    chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());
    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;


    double denominator = result.value();

    std::cout << "int mu^4 exp(-H/kT) / int \mu^2 exp(-H/kT) = " << numerator / denominator << std::endl;

    return 0;
}

