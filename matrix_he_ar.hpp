#pragma once

#include <cmath>
#include <vector>

#include "ar_he_pot.hpp"
#include "ar_he_pot_derivative.hpp"
#include "ar_he_dip_buryak_fit.hpp"

void rhs( double * out, const double R, const double pR, const double theta,
          const double pTheta, const double phi, const double pPhi );
void transform_dipole( std::vector<double> &res, const double R, const double theta, const double phi );
