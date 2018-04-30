#include <cmath>

#include "diatomic_equations_3d.hpp"
#include "ar_he_pes.hpp"
#include "ar_he_pes_der.hpp"

const double mu = 6632.039; // He-Ar

void rhs_3D(double* out, double R, double pR, double theta, double pTheta, double phi, double pPhi )
// input:
//      out -- prepared array to fill the right-hand sides of ODEs
{
    out[0] = pR / mu; // dot{R}
    out[1] = pTheta * pTheta / mu / R / R / R + \
            pPhi * pPhi / mu / R / R / R / std::sin(theta) / std::sin(theta) - ar_he_pot_der(R); // dot{pR}
    out[2] = pTheta / mu / R / R; // dot{theta}
    out[3] = pPhi * pPhi * std::cos(theta) / mu / R / R / std::pow(std::sin(theta), 3); // dot{pTheta}
    out[4] = pPhi / mu / R / R / std::sin(theta) / std::sin(theta); // dot{phi}
    out[5] = 0; // dot{pPhi}
}


