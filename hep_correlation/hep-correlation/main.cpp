#include "hep/mc.hpp"

#include "ar_he_dip_buryak_fit.hpp"
#include "ar_he_pot.hpp"
#include "ar_he_pot_derivative.hpp"
#include "constants.hpp"

#include <fstream>
#include <ctime>
#include <iostream>
#include <vector>
#include <cmath>

#include "awp.hpp"
#include "basis.hpp"
#include "gear.hpp"
#include "vmblock.hpp"

const double MU = 6632.039;
const int MaxTrajectoryLength = 500;
const double Temperature = 295.0;
const double sampling_time = 200.0; // a.t.u

const double Rint_min = 3.0;
const double Rint_max = 40.0;

void transform_dipole( std::vector<double> & res, const double R, const double theta, const double phi )
{
    double dipz = ar_he_dip_buryak_fit( R );

    res[0] = dipz * std::sin(theta) * std::cos(phi);
    res[1] = dipz * std::sin(theta) * std::sin(phi);
    res[2] = dipz * std::cos(theta);
}

void rhs( double * out, const double R, const double pR,
          const double theta, const double pTheta,
          const double phi, const double pPhi )
{
    out[0] = pR / MU; // dot{R}
    out[1] = pTheta * pTheta / MU / std::pow(R, 3) + \
            pPhi * pPhi / MU / std::pow(R, 3) / std::sin(theta) / std::sin(theta) - ar_he_pot_derivative(R); // dot{pR}
    out[2] = pTheta / MU / R / R; // dot{theta}
    out[3] = pPhi * pPhi * std::cos(theta) / MU / R / R / std::pow(std::sin(theta), 3); // dot{pTheta}
    out[4] = pPhi / MU / R / R / std::pow(sin(theta), 2); // dot{phi}
    out[5] = 0; // dot{pPhi}
}

void syst (REAL t, REAL *y, REAL *f)
{
    // параметр t мы не используем, от него у нас ничего не зависит
    (void)(t); // avoid unused parameter warning
    double *out = new double[6];

    // вызываем функцию, вычисляющую правые части
    rhs( out, y[0], y[1], y[2], y[3], y[4], y[5] );

    // засовываем в нужном порядке производные
    f[0] = out[0]; // \dot{R}
    f[1] = out[1]; // \dot{p_R}
    f[2] = out[2]; // \dot{\theta}
    f[3] = out[3]; // \dot{p_\theta}
    f[4] = out[4];
    f[5] = out[5];

    delete [] out;
}

double integrand_(hep::mc_point<double> const& x, double Temperature )
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

        double dip = ar_he_dip_buryak_fit(R);
        double res = dip * dip * std::exp(-h * constants::HTOJ / (constants::BOLTZCONST * Temperature));
        return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
    }

    else
        return 0.0;
}

// x -- точка 6-мерного фазового пространства
// ncorr -- номер значения корреляции
double correlation_integrand_(hep::mc_point<double> const& x, const int ncorr )
{
    double R_new = x.point()[0];
    double pR_new = x.point()[1];
    double theta_new = x.point()[2];
    double pTheta_new = x.point()[3];
    double phi_new = x.point()[4];
    double pPhi_new = x.point()[5];

    double R = std::tan( M_PI / 2 * R_new );
    double Theta = M_PI * theta_new;
    double Phi = 2 * M_PI * phi_new;
    double pR = std::tan(M_PI * (pR_new - 0.5));
    double pTheta = std::tan(M_PI * (pTheta_new - 0.5));
    double pPhi = std::tan(M_PI * (pPhi_new - 0.5));

    double potential_value = ar_he_pot(R);
    double h = std::pow(pR, 2) / (2 * MU) + std::pow(pTheta, 2) / (2.0 * MU * std::pow(R, 2)) + \
               std::pow(pPhi, 2.0) / (2.0 * MU * std::pow(R, 2.0) * std::pow(std::sin(Theta), 2.0)) + \
               potential_value;

    if ( R > Rint_max || R < Rint_min )
        return 0.0;

    std::vector<REAL> y0 = {R, pR, Theta, pTheta, Phi, pPhi};
    //std::cout << "R = " << R << "; pR = " << pR << std::endl;

    if ( h > 0 )
    {
        // jacobians
        double jacR = M_PI / 2.0 * (1.0 + std::pow(R, 2.0));
        double jacTheta = M_PI;
        double jacPhi = 2.0 * M_PI;
        double jacpR = M_PI * (1.0 + std::pow(pR, 2.0));
        double jacpTheta = M_PI * (1.0 + std::pow(pTheta, 2.0));
        double jacpPhi = M_PI * (1.0 + std::pow(pPhi, 2.0));

        double jacs = jacR * jacpR * jacTheta * jacPhi * jacpTheta * jacpPhi;

        /*
        // -----------------------------------------------
        REAL epsabs = 1E-13; // absolute error bound
        REAL epsrel = 1E-13; // relative error bound

        REAL t0 =  0;
        REAL h = 0.1; // initial, final step size
        REAL xend = sampling_time; // right edge of integration interval

        long fmax = 1e9; // maximum number of calls of right side in gear4()
        long aufrufe; // actual number of function calls
        int fehler;  // error code from gear4()
        // -----------------------------------------------

        std::vector<double> initial_dipole( 3 );
        transform_dipole(initial_dipole, R, Theta, Phi);

        for ( int counter = 0; counter < ncorr; counter++ )
        {
            fehler = gear4( &t0, xend, 6, syst, &y0[0], epsabs, epsrel, &h, fmax, &aufrufe );
            if ( fehler != 0 )
            {
                std::cout << "Gear4: error n = " << 10 + fehler << std::endl;
                break;
            }

            xend = sampling_time * (counter + 2);
            aufrufe = 0;
        }


        std::vector<double> final_dipole( 3 );
        transform_dipole( final_dipole, y0[0], y0[2], y0[4] );

        double prod = initial_dipole[0] * final_dipole[0] + \
                      initial_dipole[1] * final_dipole[1] + \
                      initial_dipole[2] * final_dipole[2];
        */

        double prod = ar_he_dip_buryak_fit(R) * ar_he_dip_buryak_fit(R);

        return std::exp(-h * constants::HTOJ / (constants::BOLTZCONST * Temperature)) * prod * jacs;
    }
    else
        return 0.0;

}

double denominator(hep::mc_point<double> const& x )
{
    double R_new = x.point()[0];
    double pR_new = x.point()[1];
    double theta_new = x.point()[2];
    double pTheta_new = x.point()[3];
    double Phi_new = x.point()[4];
    double pPhi_new = x.point()[5];

    double R = std::tan( M_PI / 2 * R_new );
    double Theta = M_PI * theta_new;
    double Phi = 2 * M_PI * Phi_new;

    double pR = std::tan(M_PI * (pR_new - 0.5));
    double pTheta = std::tan(M_PI * (pTheta_new - 0.5));
    double pPhi = std::tan(M_PI * (pPhi_new - 0.5));

    if ( (R > Rint_max) || (R < Rint_min) )
        return 0.0;

    double potential_value = ar_he_pot(R);
    double h = std::pow(pR, 2) / (2 * MU) + std::pow(pTheta, 2) / (2.0 * MU * std::pow(R, 2)) + \
               std::pow(pPhi, 2.0) / (2.0 * MU * std::pow(R, 2.0) * std::pow(std::sin(Theta), 2.0)) + \
               potential_value;

    if ( h > 0 )
    {
        double jacR = M_PI / 2.0 * (1.0 + std::pow(R, 2.0));
        double jacTheta = M_PI;
        double jacPhi = 2.0 * M_PI;
        double jacpR = M_PI * (1.0 + std::pow(pR, 2.0));
        double jacpTheta = M_PI * (1.0 + std::pow(pTheta, 2.0));
        double jacpPhi = M_PI * (1.0 + std::pow(pPhi, 2.0));

        double res = std::exp(-h * constants::HTOJ / (constants::BOLTZCONST * Temperature));
        return res * jacR * jacTheta * jacPhi * jacpR * jacpTheta * jacpPhi;
    }

    else
        return 0.0;
}

int main()
{
    std::clock_t start = std::clock();

    // set the verbose vegas callback function
    hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

    auto integrand = bind( correlation_integrand_, std::placeholders::_1, 0 );

    auto test = bind( integrand_, std::placeholders::_1, 295.0 );

    auto results = hep::vegas(
        hep::make_integrand<double>(integrand, 6),
        std::vector<std::size_t>(5, 1e4)
    );

    auto result = hep::cumulative_result0(results.begin() + 1, results.end());
    double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;

    std::cout << "Time elapsed: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << " s" << std::endl;

    double numerator = result.value();

    results = hep::vegas(
         hep::make_integrand<double>(denominator, 6),
         std::vector<std::size_t>(5, 1e4)
                );

    result = hep::cumulative_result0(results.begin() + 1, results.end());

    chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());
    std::cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls()
              << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << std::endl;

    double denominator = result.value();

    double res = numerator / denominator;
    std::cout << "ratio = " << res << std::endl;

    return 0;
}

