#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

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
const double Rint_min = 3.0; // a0
const double Rint_max = 25.0; // a0
// ------------------------------------------------------

class SplinedData
{
private:
    static std::vector<double> j_values;
    static std::vector<double> xmax_values;

public:
    SplinedData( std::string const & filename )
    {
        std::ifstream spline_file(filename);

        if ( !spline_file )
            throw std::invalid_argument("Can't find splined file.");

        const int MAX_LINE = 100;
        char line[MAX_LINE];

        while ( spline_file.getline(line, MAX_LINE) )
        {
            std::stringstream ss;
            ss.str(line);

            double j_value;
            ss >> j_value;

            double xmax_value;
            ss >> xmax_value;

            j_values.push_back( j_value );
            xmax_values.push_back( xmax_value );
        }

        spline_file.close();
    }

    void show( )
    {
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Spline file is successfully read \n";
        std::cout << "Min J value: " << j_values[0] << "; max J value: " << j_values.back() << std::endl;
        std::cout << "Min Xmax value: " << xmax_values[0] << "; max Xmax value: " << xmax_values.back() << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
    }

    double xmax( double j ) const
    {
        // нет максимума (либо прыгнули выше перегиба эффективного потенциала или слишком маленький угловой момент
        if  ( j < j_values[0] || j > j_values.back() )
        {
            //std::cout << "j less than" << j_values[0] << " or j more than " << j_values.back() << std::endl;
            return -1;
        }

        // нижняя граница для полученного значения j в сортированном массиве j_values
        std::vector<double>::iterator up_it = std::upper_bound( j_values.begin(), j_values.end(), j );
        size_t up = up_it - j_values.begin();

        // определяем верхнюю границу для полученного значения j в массиве j_values
        size_t low = up - 1;

        //std::cout << "input: " << j << "; low value: " << j_values[low] << "; low xmax: " << xmax_values[low] << std::endl;
        //std::cout << "input: " << j << "; up value: " << j_values[up] << "; up xmax: " << xmax_values[up] << std::endl;

        // если значение в списке
        if ( j == j_values[low] )
            return xmax_values[low];

        // иначе выдаем линейную аппроксимацию между найденными границами
        double k = (xmax_values[up] - xmax_values[low]) / (j_values[up] - j_values[low]);
        //std::cout << "k: " << k << std::endl;

        return xmax_values[up] - k * (j_values[up] - j);
    }
};
// инициализация статических векторов
std::vector<double> SplinedData::j_values;
std::vector<double> SplinedData::xmax_values;

double integrand_(hep::mc_point<double> const& x, double Temperature, SplinedData const & sd, bool numerator )
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

    double j = std::sqrt( pTheta*pTheta + pPhi * pPhi / std::sin(Theta) / std::sin(Theta) );
    double xmax = sd.xmax( j ); // положение максимума
    double eff_pot = ar_he_pot(xmax) + std::pow(j, 2.0) / (2.0 * MU * std::pow(xmax, 2));

    // квазисвязанное состояние: левее локального максимума и энергия меньше максимума
    if ( R < xmax && h < eff_pot )
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
        // в числителе стоит интеграл \mu^2 \exp(-H/kT)
        if ( numerator )
        {
            double dip = ar_he_dip_buryak_fit(R);
            double res = dip * dip * std::exp(-h * HTOJ / (BOLTZCONST * Temperature));
            return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
        }
        // в знаменателе стоит интеграл \exp(-H/kT)
        else
        {
            double res = std::exp(-h * HTOJ / (BOLTZCONST * Temperature));
            return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
        }
    }

    else
        return 0.0;
}

int main()
{
    SplinedData sd("../j_xmax_values.txt");
    sd.show();

    const double Temperature = 295.0;
    std::cout << "Current temperature = " << Temperature << std::endl;

    // set the verbose vegas callback function
    hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

    auto num_integrand = bind( integrand_, std::placeholders::_1, Temperature, sd, true );
    auto denom_integrand = bind( integrand_, std::placeholders::_1, Temperature, sd, false );

    auto results = hep::vegas(
       hep::make_integrand<double>(num_integrand, 6),
       std::vector<std::size_t>(10, 2e6)
     );

    auto result = hep::cumulative_result0(results.begin() + 1, results.end());
    double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());
    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;

    double numerator = result.value();

    results = hep::vegas(
         hep::make_integrand<double>(denom_integrand, 6),
         std::vector<std::size_t>(10, 2e6)
                );

    result = hep::cumulative_result0(results.begin() + 1, results.end());

    chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());
    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;


    double denominator = result.value();

    std::cout << "int mu^2 exp(-H/kT) / int exp(-H/kT) = " << numerator / denominator << std::endl;

    double Volume = 4.0 / 3.0 * M_PI * std::pow(Rint_max, 3.0) * std::pow(ALU, 3);

    std::cout << " int/int x Volume = " << numerator / denominator * Volume << std::endl;

    return 0;
}
