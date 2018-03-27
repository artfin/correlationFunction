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

int Trajectory::report_trajectory_status( void )
{
    if ( cut_trajectory )
    {
        std::cout << "(slave): sending cutting trajectory signal" << std::endl;
        MPI_Send( &cut_trajectory, 1, MPI_INT, 0, tags::TRAJECTORY_CUT_TAG, MPI_COMM_WORLD );
    }
    else
    {
        MPI_Send( &cut_trajectory, 1, MPI_INT, 0, tags::TRAJECTORY_FINISHED_TAG, MPI_COMM_WORLD );
    }

    return cut_trajectory;
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

    // первый шаг -- половинный по времени
    //REAL xend = parameters.sampling_time / 2.0;

    long fmax = 1e9; // maximum number of calls of right side in gear4()
    long aufrufe; // actual number of function calls
    int fehler;  // error code from gear4()

    std::vector<double> temp( 3 );

    int counter = 1;
    while( y0[0] < parameters.RDIST )
    {
        if ( counter == parameters.MaxTrajectoryLength / 2 )
        {
            std::cout << "Trajectory cut!" << std::endl;
            cut_trajectory = 1;
            break;
        }

        fehler = gear4( &t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe );
        if ( fehler != 0 )
        {
            std::cout << "Gear4: error n = " << 10 + fehler << std::endl;
            break;
        }

        if ( counter % 1000 == 0 )
        {
            std::cout << "..." << std::endl;
        }

        //std::vector<double> coords{ t0, y0[0], y0[1], y0[2], y0[3] };
        //trajectory.push_back( coords );

        transform_dipole( temp, y0[0], y0[2], y0[4] );
        dipx.push_back( temp[0] );
        dipy.push_back( temp[1] );
        dipz.push_back( temp[2] );

        // точки вида (i + 1/2) * sampling_time
        xend = parameters.sampling_time * counter + parameters.sampling_time / 2.0;
        aufrufe = 0;

        counter++;
    }
}

