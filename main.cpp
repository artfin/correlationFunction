//#include <mpi.h>
#include <mpi/mpi.h>

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

const int CORRELATION_FUNCTION_LENGTH = 3000;

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

void save_correlation_function_final( std::vector<double> const & v, const std::string filename, const double sampling_time )
{
    double time = 0.0;

    std::ofstream outFile( filename );
    for ( size_t k = 0; k < v.size(); k++, time += sampling_time )
        outFile << time << " " << v[k] << std::endl;

    outFile.close();
}

void save_correlation_function_debug( std::vector<double> const & v, const std::string filename, const double constant )
{
    std::ofstream outFile( filename );
    for ( size_t k = 0; k < v.size(); k++ )
        outFile << v[k] * constant << std::endl;

    outFile.close();
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

    // sending first trajectory
    for ( int i = 1; i < world_size; i++ )
    {
        p = MHA_generate_point( p, getOutOf, density );
        addPoint( histograms, p );

        MPI_Send( &p[0], DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent point " << endl;

        MPI_Send( &sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent number of trajectory" << endl;

        sent++;
    }

    std::vector<double> correlation_total(CORRELATION_FUNCTION_LENGTH);
    std::vector<double> correlation_package(CORRELATION_FUNCTION_LENGTH);

    const int toSave = 1000; // раз в какое количество точек сохранять файл
    std::string filenameCorrelationFunction = "../output/equilibrium_correlation.txt";

    double Volume = 4.0 / 3.0 * M_PI * pow( Rint_max * constants::ALU, 3);

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
            std::cout << "(master) Received cutting trajectory tag!" << std::endl;
            if ( sent <= parameters.NPOINTS )
            {
                p = MHA_generate_point( p, getOutOf, density );
                addPoint( histograms, p );

                MPI_Send( &p[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
                MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
                //cout << "(master) Sent new point" << endl;

                continue;
            }
        }

        std::fill( correlation_package.begin(), correlation_package.end(), 0.0 );

        MPI_Recv( &correlation_package[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

        for ( size_t i = 0; i < CORRELATION_FUNCTION_LENGTH; i++ )
        {
            correlation_total[i] += correlation_package[i];
        }

        received++;
        //cout << "(master) after all MPI_Recv; received = " << received << endl;

        if ( (received + 1) % toSave == 0 )
        {
            std::cout << "(debug) Saving correlation file to " << filenameCorrelationFunction << std::endl;
            save_correlation_function_debug( correlation_total, filenameCorrelationFunction, Volume / received );
        }

        if ( received == parameters.NPOINTS )
        {
            for ( size_t k = 0; k < correlation_total.size(); k++ )
            {
                correlation_total[k] *= Volume / received;
            }

            std::cout << "(final) Saving correlation file to " << filenameCorrelationFunction << std::endl;
            save_correlation_function_final( correlation_total, filenameCorrelationFunction,
                                             parameters.sampling_time * constants::ATU );
            is_finished = true;
        }

        if ( sent < parameters.NPOINTS )
        {
            p = MHA_generate_point( p, getOutOf, density );
            addPoint( histograms, p );

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
        // нет встроенного типа MPI_BOOL. )
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

        // начинаем расчет траектории
        trajectory.run_trajectory( syst );

        // собираем статус траектории. метод отправляет статус траектории мастеру
        // если траектория обрублена, то обнуляем буферы хранения диполя с траектории
        cut_trajectory = trajectory.report_trajectory_status( );
        if ( cut_trajectory )
        {
            trajectory.dump_dipoles();
            continue;
        }

        // если мы прошли предыдущий блок кода, значит траектория имеет допустимую длину.
        // копируем компоненты дипольного момента вдоль траектории в виде структуры данных vector<double>
        dipx_forward = trajectory.get_dipx();
        dipy_forward = trajectory.get_dipy();
        dipz_forward = trajectory.get_dipz();

        // после копирования освобождаем эти вектора внутри объекта trajectory
        trajectory.dump_dipoles( );

        std::cout << "(" << world_rank << ") Processing " << trajectory.get_trajectory_counter()
                  << " trajectory. npoints = " << dipz_forward.size() << "; time = "
                  << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

        calculate_physical_correlation(correlation_function, dipx_forward, dipy_forward, dipz_forward);

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
