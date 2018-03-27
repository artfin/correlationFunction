#include "matrix_he_ar.hpp"

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
    const double MU = 6632.039;

    out[0] = pR / MU; // dot{R}
    out[1] = pTheta * pTheta / MU / std::pow(R, 3) + \
            pPhi * pPhi / MU / std::pow(R, 3) / std::sin(theta) / std::sin(theta) - ar_he_pot_derivative(R); // dot{pR}
    out[2] = pTheta / MU / R / R; // dot{theta}
    out[3] = pPhi * pPhi * std::cos(theta) / MU / R / R / std::pow(std::sin(theta), 3); // dot{pTheta}
    out[4] = pPhi / MU / R / R / std::pow(sin(theta), 2); // dot{phi}
    out[5] = 0; // dot{pPhi}
}
