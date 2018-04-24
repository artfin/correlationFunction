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

const int DIM = 6;
const int CORRELATION_FUNCTION_LENGTH = 3000;

const double Rint_max = 40.0;

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

void save_correlation_function( std::vector<double> const & v, const std::string filename,
                                const double sampling_time, const double normalizing_constant )
{
    std::ofstream outFile(filename);

    double time = 0.0;
    for ( size_t k = 0; k < v.size(); ++k, time += sampling_time )
        outFile << time << " " << v[k] * normalizing_constant << std::endl;

    outFile.close();
}

void read_initial_conditions( std::vector<std::vector<double>> & contents, const std::string filename,
                              int lines_to_read = -1 )
{
    std::ifstream inFile(filename);
    const int MAXLINE = 100;
    char buf[MAXLINE];

    std::stringstream ss;
    std::vector<double> temp;
    double value;

    int lines_counter = 0;
    while ( inFile.getline(buf, MAXLINE) )
    {
        ss << buf;
        for ( int k = 0; k < DIM; ++k )
        {
            ss >> value;
            temp.push_back(value);
        }

        if ( lines_counter == lines_to_read )
            break;

        contents.push_back(temp);
        temp.clear();
        ss.clear();
        ss.str("");

        ++lines_counter;
    }

    inFile.close();
}

void merge( std::vector<double> & full, std::vector<double> & backward, std::vector<double> & forward )
{
    full.resize(backward.size() + forward.size());
    std::reverse_copy(backward.begin(), backward.end(), full.begin());
    std::copy(forward.begin(), forward.end(), full.begin() + backward.size());
}

void master_code( int world_size )
{
    MPI_Status status;
    int source;

    Parameters parameters;
    FileReader fileReader( "../parameters.in", &parameters );

    std::vector<std::vector<double>> samples;
    const int samples_to_read = 1000;
    read_initial_conditions(samples, "../initial_points_transformed_space_1e5.txt", samples_to_read);
    std::cout << "Samples len: " << samples.size() << std::endl;

    const int toSave = 1000; // размер сохраняемого блока
    assert( samples.size() % toSave == 0 );

    int samples_counter = 0;

    int sent = 0;
    int received = 0;

    // status of calculation
    bool is_finished = false;

    std::vector<double> sample;

    // sending first trajectory
    for ( int i = 1; i < world_size; i++ )
    {
        sample = samples[samples_counter];
        ++samples_counter;

        MPI_Send( &sample[0], DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent point " << endl;

        MPI_Send( &sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent number of trajectory" << endl;

        ++sent;
    }

    std::vector<double> correlation_block(CORRELATION_FUNCTION_LENGTH);
    std::vector<double> correlation_total(CORRELATION_FUNCTION_LENGTH);
    std::vector<double> correlation_package(CORRELATION_FUNCTION_LENGTH);

    int filesCounter = 0; // номер файла
    std::string filenameCorrelationFunction = "../output/eqcorr";

    double Volume = 4.0 / 3.0 * M_PI * pow( Rint_max * constants::ALU, 3);

    int total_points = 0;
    int traj_len = 0;

    while( true )
    {
        //std::cout << "(master) cycle beginning" << std::endl;

        if ( is_finished )
        {
            for ( int i = 1; i < world_size; i++ )
                MPI_Send( &is_finished, 1, MPI_INT, i, tags::EXIT_TAG, MPI_COMM_WORLD );

            break;
        }

        std::fill( correlation_package.begin(), correlation_package.end(), 0.0 );
        MPI_Recv( &correlation_package[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        source = status.MPI_SOURCE;

        MPI_Recv( &traj_len, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        total_points += traj_len;

        for ( size_t i = 0; i < CORRELATION_FUNCTION_LENGTH; i++ )
        {
            correlation_total[i] += correlation_package[i];
            correlation_block[i] += correlation_package[i];
        }
        ++received;


        if ( (received % toSave == 0)  && (received != 0) )
        {
            std::string filename = filenameCorrelationFunction + std::to_string(filesCounter) + ".txt";
            ++filesCounter;

            std::cout << "(debug) Saving correlation file to " << filename << std::endl;
            //save_correlation_function( correlation_block, filename,
            //                          parameters.sampling_time * constants::ATU, Volume / toSave );
            save_correlation_function( correlation_block, filename,
                                       parameters.sampling_time * constants::ATU, 1.0 / toSave );

            // обнуляем блок
            std::fill( correlation_block.begin(), correlation_block.end(), 0.0 );
        }

        if ( received == (int) samples.size() )
        {
            std::string filename = filenameCorrelationFunction + "_final.txt";
            std::cout << "(final) Saving final correlation to " << filename << std::endl;
            //save_correlation_function( correlation_total, filename,
            //                           parameters.sampling_time * constants::ATU, Volume / received );
            save_correlation_function( correlation_total, filename,
                                       parameters.sampling_time * constants::ATU, 1.0 / received );


            std::cout << "(final) Total points (on trajectories): " << total_points << std::endl;
            std::cout << "(final) Exiting..." << std::endl;
            is_finished = true;
        }

        if ( sent < (int) samples.size() )
        {
            sample = samples[samples_counter];
            ++samples_counter;

            MPI_Send( &sample[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
            MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );

            ++sent;
        }
    }
}

// weight = R^2 \sin \Theta
void calculate_physical_correlation( std::vector<double> & physical_correlation,
                                     std::vector<double> & dipx,
                                     std::vector<double> & dipy,
                                     std::vector<double> & dipz, const double weight )
{
    size_t size = dipx.size();
    assert( dipy.size() == size );
    assert( dipz.size() == size );

    int max_size = std::min(CORRELATION_FUNCTION_LENGTH, (int) size);

    // подготавливает вектор нужного размера, заполняет его элементы нулями
    physical_correlation.resize( CORRELATION_FUNCTION_LENGTH );


    // E = mu(0)*mu(t) / mu(0)*mu(0)
    // без трюка
    double dip_initial = dipx[0] * dipx[0] + dipy[0] * dipy[0] + dipz[0] * dipz[0];
    for ( int n = 0; n < max_size; ++n )
    {
        physical_correlation[n] = ( dipx[0] * dipx[n] + \
                                    dipy[0] * dipy[n] + \
                                    dipz[0] * dipz[n] ) / dip_initial * weight;
    }

    // с трюком
    /*
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
    */
}

void slave_code( int world_rank )
{
    Parameters parameters;
    FileReader fileReader( "../parameters.in", &parameters );

    bool exit_status = false;
    Trajectory trajectory( parameters );

    //std::vector<double> dipx_backward, dipy_backward, dipz_backward;
    std::vector<double> dipx_forward, dipy_forward, dipz_forward;
    //std::vector<double> dipx, dipy, dipz;

    std::vector<double> correlation_function;

    int trajectory_len;


    while ( true )
    {
        std::clock_t start = clock();

        // реализуем MPI получение начальных условий и запись их в private элементы объекта класса Trajectory
        // возвращаем статус полученного сообщения. Если статус полученного сообщения таков, что больше траекторий
        // обсчитывать не надо, то exit_status = true и, соответственно, текущий процесс выходит из бесконечного цикла
        // и прекращает свою работу.
        exit_status = trajectory.receive_initial_conditions( );

        double weight = trajectory.get_weight();

        if ( exit_status )
        {
            std::cout << "(" << world_rank << ") exit message is received!" << std::endl;
            break;
        }

        // начинаем расчет траектории
        trajectory.run_trajectory( syst );

        // если мы прошли предыдущий блок кода, значит траектория имеет допустимую длину.
        // копируем компоненты дипольного момента вдоль траектории в виде структуры данных vector<double>
        dipx_forward = trajectory.get_dipx();
        dipy_forward = trajectory.get_dipy();
        dipz_forward = trajectory.get_dipz();

        // после копирования освобождаем эти вектора внутри объекта trajectory
        trajectory.dump_dipoles( );

        /*
        trajectory.reverse_initial_conditions();

        trajectory.run_trajectory( syst );

        dipx_backward = trajectory.get_dipx();
        dipy_backward = trajectory.get_dipy();
        dipz_backward = trajectory.get_dipz();

        trajectory.dump_dipoles();

        merge(dipx, dipx_backward, dipx_forward);
        merge(dipy, dipy_backward, dipy_forward);
        merge(dipz, dipz_backward, dipz_forward);

        dipx_backward.clear(); dipx_forward.clear();
        dipy_backward.clear(); dipy_forward.clear();
        dipz_backward.clear(); dipz_forward.clear();
        */

        std::cout << "(" << world_rank << ") Processing " << trajectory.get_trajectory_counter()
                  << " trajectory. npoints = " << dipz_forward.size() << "; time = "
                  << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

        calculate_physical_correlation(correlation_function, dipx_forward, dipy_forward, dipz_forward, weight);

        // Отправляем собранный массив корреляций мастер-процессу
        MPI_Send( &correlation_function[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        trajectory_len = (int) dipz_forward.size();
        MPI_Send( &trajectory_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
        correlation_function.clear();

        dipx_forward.clear();
        dipy_forward.clear();
        dipz_forward.clear();
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
