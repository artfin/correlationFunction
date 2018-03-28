#include <mpi.h>

#include <iostream>
#include <ctime>
#include <chrono>
#include <cassert>
#include <random>
#include <functional>
#include <algorithm>

#include "constants.hpp"
#include "filereader.hpp"
#include "parameters.hpp"
#include "ar_he_pot.hpp"

#include "trajectory.hpp"

unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 mcmc_generator{ seed };

const int CORRELATION_FUNCTION_LENGTH = 100;

const double step = 1.5;
const double DIM = 6;

const double Racc_max = 40.0;
const double Racc_min = 4.0;

const double MU = 6632.039;

void sample(std::vector<double> const & prev, std::vector<double> & next)
{
    std::normal_distribution<double> distribution( 0.0, step );
    for ( size_t i = 0; i < next.size(); ++i )
         next[i] = prev[i] + distribution( mcmc_generator );
}

// порядок перменных: R, pR, theta, pT, phi, pPhi
double density_( std::vector<double> & x, const double Temperature )
{
    double R = x[0];
    double pR = x[1];
    double theta = x[2];
    double pTheta = x[3];
    // double phi = x[4];
    double pPhi = x[5];

    double H = std::pow(pR, 2) / (2.0 * MU) + \
               std::pow(pTheta, 2) / (2.0 * MU * R * R) + \
               std::pow(pPhi, 2) / (2.0 * MU * R * R * sin(theta) * sin(theta)) + \
               ar_he_pot(R);

    if ( (H > 0) && (R > Racc_min) && (R < Racc_max) )
        return exp(- H * constants::HTOJ / (constants::BOLTZCONST * Temperature) );
    else
        return 0.0;
}

std::vector<double> MHA_burnin( int burnin, std::vector<double> x, std::function<double(std::vector<double>&)> density )
{
    std::vector<double> xcand(DIM);
    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    for (int i = 0; i < burnin; ++i)
    {
        sample(x, xcand);
        double alpha = density(xcand) / density(x);

        if ( unidistr(mcmc_generator) < alpha)
        {
            x = xcand;
        }
    }

    return x;
}

// getOutOf -- то количество точек, раз в которое мы выбираем точку,
// чтобы на ней посчитать интегранд
std::vector<double> MHA_generate_point(std::vector<double> x, int getOutOf,
                                       std::function<double(std::vector<double>&)> density)
{
    std::vector<double> xcand(DIM);
    int counter = 1;

    std::uniform_real_distribution<double> unidistr(0.0, 1.0);

    while( true)
    {
        sample(x, xcand);
        double alpha = density(xcand) / density(x);

        if ( unidistr(mcmc_generator) < std::min(alpha, 1.0) )
        {
            x = xcand;
        }

        if ( (counter > getOutOf) && (xcand[0] > Racc_min) && (xcand[0] < Racc_max) )
        {
            return xcand;
        }

        ++counter;
    }
}

// "интерфейсная" функция. она передается методу GEAR, который
// осуществляет решение системы дифуров, правая часть которой, задается этой
// функцией. Внутри мы осуществляем вызов функции, которая лежит в файле
// matrix_he_ar.cpp, там осуществлен матричный расчет правой части дифуров.
// Эта функция сделана для удобства, устроена таким образом, чтобы было
// совместимо с gear4.
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

void master_code( int world_size )
{
    MPI_Status status;
    int source;

    Parameters parameters;
    FileReader fileReader( "../parameters.in", &parameters );

    int sent = 0;
    int received = 0;

    // status of calculation
    bool is_finished = false;

    std::function<double(std::vector<double>&)> density = bind( density_, std::placeholders::_1, parameters.Temperature );

    const int burnin = 3e4;
    std::vector<double> initial_after_burnin = MHA_burnin( burnin, parameters.initial_point, density );

    const int getOutOf = 3;
    std::vector<double> p = MHA_generate_point( initial_after_burnin, getOutOf, density );

    // sending first trajectory
    for ( int i = 1; i < world_size; i++ )
    {
        p = MHA_generate_point( p, getOutOf, density );

        MPI_Send( &p[0], DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent point " << endl;

        MPI_Send( &sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent number of trajectory" << endl;

        sent++;
    }

    std::clock_t start = clock();

    std::vector<double> correlation_total(CORRELATION_FUNCTION_LENGTH);
    std::vector<double> correlation_package(CORRELATION_FUNCTION_LENGTH);

    while( true )
    {
        if ( is_finished )
        {
            for ( int i = 1; i < world_size; i++ )
                MPI_Send( &is_finished, 1, MPI_INT, i, tags::EXIT_TAG, MPI_COMM_WORLD );

            break;
        }

        int msg;
        MPI_Recv( &msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        source = status.MPI_SOURCE;
        if ( status.MPI_TAG == tags::TRAJECTORY_CUT_TAG )
        {
            //cout << "(master) Received cutting trajectory tag!" << endl;
            if ( sent <= parameters.NPOINTS )
            {
                p = MHA_generate_point( p, getOutOf, density );

                MPI_Send( &p[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
                MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
                //cout << "(master) Sent new point" << endl;

                continue;
            }
        }

        /*
        if ( status.MPI_TAG == tags::TRAJECTORY_FINISHED_TAG )
            cout << "(master) Trajectory is not cut!" << endl;
        */

        std::fill( correlation_package.begin(), correlation_package.end(), 0.0 );

        MPI_Recv( &correlation_package[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

        double Volume = 4.0 / 3.0 * M_PI * pow(parameters.RDIST * constants::ALU, 3);
        double correlationFunctionConstant = Volume; // * constants::ADIPMOMU * constants::ADIPMOMU \
                                / (4 * M_PI * constants::EPSILON0);

        if ( received >= 1 )
        {
            for ( size_t i = 0; i < CORRELATION_FUNCTION_LENGTH; i++ )
            {
                // размерность корреляции дипольного момента -- квадрат диполя
                correlation_total[i] += correlation_package[i] * correlationFunctionConstant;
                correlation_total[i] *= 0.5;
            }
        }
        else
        {
            for ( size_t i = 0; i < CORRELATION_FUNCTION_LENGTH; i++ )
                correlation_total[i] += correlation_package[i] * correlationFunctionConstant;
        }

        received++;
        //cout << "(master) after all MPI_Recv; received = " << received << endl;

        if ( received == parameters.NPOINTS )
        {
            std::ofstream file( "../output/equilibrium_correlation.txt" );
            for ( size_t k = 0; k < CORRELATION_FUNCTION_LENGTH; k++ )
                file << correlation_total[k] << std::endl;
            file.close();

            is_finished = true;
        }

        if ( sent < parameters.NPOINTS )
        {
            p = MHA_generate_point( p, getOutOf, density );

            MPI_Send( &p[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
            MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );

            sent++;
        }
    }
}

void calculate_physical_correlation( std::vector<double> & physical_correlation,
                                     std::vector<double> & dipx,
                                     std::vector<double> & dipy,
                                     std::vector<double> & dipz )
{
    size_t size = dipx.size();
    assert( dipy.size() == size );
    assert( dipz.size() == size );

    int max_size = std::min(CORRELATION_FUNCTION_LENGTH, (int) size);

    // подготавливает вектор нужного размера, заполняет его элементы нулями
    physical_correlation.resize( CORRELATION_FUNCTION_LENGTH );

    double res = 0;
    for ( int n = 0; n < max_size; n++ )
    {
        res = 0;
        for ( size_t i = 0, j = n; j < size; i++, j++ )
        {
            res += dipx[i] * dipx[j];
            res += dipy[i] * dipy[j];
            res += dipz[i] * dipz[j];
        }

        physical_correlation[n] = res / (size - n);
    }
}

void slave_code( int world_rank )
{
    Parameters parameters;
    FileReader fileReader( "../parameters.in", &parameters );

    int cut_trajectory = 0;
    bool exit_status = false;
    Trajectory trajectory( parameters );

    std::vector<double> dipx_forward, dipy_forward, dipz_forward;
    std::vector<double> correlation_function;

    while ( true )
    {
        // переменная cut_trajectory играет роль переменной типа bool
        // ( это сделано для того, чтобы переменная могла быть переслана при помощи MPI_Send,
        // там нет встроенного типа MPI_BOOL. )
        //
        // если траектория оказывается обрублена (длиннее чем parameters.MaxTrajectoryLength точек),
        // то в классе Trajectory статус траектории становится cut и при помощи метода
        // .report_trajectory_status() получаем статус текущей траектории.
        // В начале каждой итерации мы поэтому должны занулить этот статус здесь и внутри объекта класса Trajectory
        cut_trajectory = 0;
        trajectory.set_cut_trajectory( 0 );

        std::clock_t start = clock();

        // реализуем MPI получение начальных условий и запись их в private элементы объекта класса Trajectory
        // возвращаем статус полученного сообщения. Если статус полученного сообщения таков, что больше траекторий
        // обсчитывать не надо, то exit_status = true и, соответственно, текущий процесс выходит из бесконечного цикла
        // и прекращает свою работу.
        exit_status = trajectory.receive_initial_conditions( );
        if ( exit_status )
        {
            std::cout << "(" << world_rank << ") exit message is received!" << std::endl;
            break;
        }

        // начинаем траекторию из полученного начального положения в фазовом пространстве
        trajectory.run_trajectory( syst );

        // если мы прошли предыдущий блок кода, значит траектория имеет допустимую длину.
        // копируем компоненты дипольного момента вдоль траектории в виде структуры данных vector<double>
        dipx_forward = trajectory.get_dipx();
        dipy_forward = trajectory.get_dipy();
        dipz_forward = trajectory.get_dipz();
        // после копирования освобождаем эти вектора внутри объекта trajectory
        trajectory.dump_dipoles( );

        cut_trajectory = trajectory.report_trajectory_status( );
        if ( cut_trajectory )
        {
            trajectory.dump_dipoles();
            continue;
        }

        std::cout << "(" << world_rank << ") Processing " << trajectory.get_trajectory_counter()
                  << " trajectory. npoints = " << dipz_forward.size() << "; time = "
                  << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

        calculate_physical_correlation(correlation_function, dipx_forward, dipy_forward, dipz_forward);
        std::cout << "Mean dipole on trajectory: " << correlation_function[0] << std::endl;

        dipx_forward.clear();
        dipy_forward.clear();
        dipz_forward.clear();

        // Отправляем собранный массив корреляций мастер-процессу
        MPI_Send( &correlation_function[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        correlation_function.clear();
    }
}

int main( int argc, char * argv[] )
{
    //Initialize the MPI environment
    MPI_Init( &argc, &argv );

    //getting id of the current process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //getting number of running processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if ( world_rank == 0 )
    {
        std::clock_t start = clock();

        master_code( world_size );

        std::cout << "Time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;
    }
    else
    {
        slave_code( world_rank );
    }

    MPI_Finalize();

    return 0;
}
