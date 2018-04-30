#include <math.h>

#include "diatomic_equations_2d.hpp"
#include "ar_he_pes.hpp"
#include "ar_he_pes_der.hpp"

const double mu = 6632.039; // He-Ar

void rhs_2D(double* out, double R, double pR, double psi, double ppsi)
// input:
//      out -- prepared array to fill the right-hand sides of ODEs
{
    out[0] = pR/mu; // /dot(R) = dH/dpR
    out[1] = ppsi*ppsi/ mu / R /R /R - ar_he_pot_der(R);  // /dot(pR) = - dH/dR = - dT/dR - dU/dR (to hartrees from cm^-1)
    out[2] = ppsi/mu/R/R; // /dot(psi) = dH/dppsi
    out[3] = 0.0;  // /dot(ppsi) = - dH/dpsi
}

