#include <math.h>

#include "diatomic_equations_2D.h"
#include "ar_he_pes.h"
#include "ar_he_pes_der.h"

//const double mu = 6631.8; //AU, He-Ar
const double mu = 6632.039;

void rhs_2D(double* out, double R, double pR, double psi, double ppsi)
// input:
//      out -- prepared array to fill the right-hand sides of ODEs
{
 //   double Jsint = J * sin(theta);

    
   // hamiltonian(derivatives, R, Theta, pR, pT, phi, theta, J);
    
    out[0] = pR/mu; // /dot(R) = dH/dpR   
    out[1] = ppsi*ppsi/ mu / R /R /R - ar_he_pot_der(R);  // /dot(pR) = - dH/dR = - dT/dR - dU/dR (to hartrees from cm^-1)
    out[2] = ppsi/mu/R/R; // /dot(psi) = dH/dppsi
    out[3] = 0.0;  // /dot(ppsi) = - dH/dpsi 
  

   // out[4] = pPhi/mu/R/R/pow(sin(theta), 2.0); // /dot(phi) = dH/dpPhi
    //out[5] = 0.0 ;  // /dot(pPhi) = - dH/dphi 
    
    

   
}
