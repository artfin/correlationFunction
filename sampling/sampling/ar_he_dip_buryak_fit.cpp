#include "ar_he_dip_buryak_fit.hpp"

double ar_he_dip_buryak_fit( const double & R )
{
    return std::exp( 0.47293181 - 0.40100613 * R - 0.10292726 * R * R ) -148.55391 * std::pow( R, -7 );
}

