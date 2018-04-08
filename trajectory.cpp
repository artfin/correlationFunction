#include "trajectory.hpp"

bool Trajectory::receive_initial_conditions( void )
{
    MPI_Status status;
    MPI_Recv( y0, parameters.DIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

    if ( status.MPI_TAG == tags::EXIT_TAG )
    {
        std::cout << "Received exit tag." << std::endl;
        return true;
    }

    MPI_Recv( &trajectory_counter, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

    memcpy( y0_copy, y0, N * sizeof(REAL) );

    return false;
}

double Trajectory::hamiltonian( const double * x )
{
    const double MU = 6632.039;

    double R = x[0];
    double pR = x[1];
    double theta = x[2];
    double pTheta = x[3];
    // double phi = x[4];
    double pPhi = x[5];

    return std::pow(pR, 2) / (2.0 * MU) + \
           std::pow(pTheta, 2) / (2.0 * MU * R * R) + \
           std::pow(pPhi, 2) / (2.0 * MU * R * R * sin(theta) * sin(theta)) + \
           ar_he_pot(R);
}


void Trajectory::save_trajectory( std::string filename )
{
    std::ofstream file( filename );

    file << std::setprecision( 8 );
    for ( size_t i = 0; i != trajectory.size(); ++i)
    {
        for ( size_t j = 0; j != trajectory[i].size(); ++j )
            file << trajectory[i][j] << " ";
        file << std::endl;
    }

    trajectory.clear();
    file.close();
}

void Trajectory::reverse_initial_conditions( void )
{
    memcpy( y0, y0_copy, N * sizeof(REAL) );

    for ( size_t i = 0; i != N; i++ )
        if ( i % 2 == 1 )
            y0[i] = -y0[i];
}

void Trajectory::dump_dipoles( void )
{
    dipx.clear();
    dipy.clear();
    dipz.clear();
}

void Trajectory::set_initial_conditions( std::vector<double>& ic )
{
    memcpy( y0, &ic[0], ic.size() * sizeof(double) );
    memcpy( y0_copy, &ic[0], ic.size() * sizeof(double) );
}

void Trajectory::show_initial_conditions( void )
{
    std::cout << "y0: ";
    for ( size_t i = 0; i != N; i++ )
        std::cout << y0[i] << " ";
    std::cout << std::endl;
}

void Trajectory::run_trajectory( dglsysfnk syst )
{
    REAL epsabs = 1E-13; // absolute error bound
    REAL epsrel = 1E-13; // relative error bound

    REAL t0 =  0;
    REAL h = 0.1; // initial, final step size
    REAL xend = 1e-10; // right edge of integration interval

    long fmax = 1e9; // maximum number of calls of right side in gear4()
    long aufrufe; // actual number of function calls
    int fehler;  // error code from gear4()

    std::vector<double> temp( 3 );

    int counter = 1;
    while( y0[0] < parameters.RDIST )
    {
        if ( counter == parameters.MaxTrajectoryLength / 2 )
            break;

        fehler = gear4( &t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe );
        if ( fehler != 0 )
        {
            std::cout << "Gear4: error n = " << 10 + fehler << std::endl;
            break;
        }

        //std::vector<double> coords{ t0, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], ar_he_dip_buryak_fit(y0[0]) };
        //trajectory.push_back( coords );

        transform_dipole( temp, y0[0], y0[2], y0[4] );
        dipx.push_back( temp[0] );
        dipy.push_back( temp[1] );
        dipz.push_back( temp[2] );

        xend = parameters.sampling_time * counter;
        aufrufe = 0;

        counter++;
    }
}

