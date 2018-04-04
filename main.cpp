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

/*
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
*/

void save_correlation_function( std::vector<double> const & v, const std::string filename,
                                const double sampling_time, const double normalizing_constant )
{
    std::ofstream outFile(filename);

    double time = 0.0;
    for ( size_t k = 0; k < v.size(); ++k, time += sampling_time )
        outFile << time << " " << v[k] * normalizing_constant << std::endl;

    outFile.close();
}

void read_initial_conditions( std::vector<std::vector<double>> & contents, const std::string filename )
{
    std::ifstream inFile(filename);
    const int MAXLINE = 100;
    char buf[MAXLINE];

    std::stringstream ss;
    std::vector<double> temp;
    double value;

    while ( inFile.getline(buf, MAXLINE) )
    {
        ss << buf;
        for ( int k = 0; k < DIM; ++k )
        {
            ss >> value;
            temp.push_back(value);
        }

        contents.push_back(temp);
        temp.clear();
        ss.clear();
        ss.str("");
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
    read_initial_conditions(samples, "../samples.txt");
    std::cout << "Samples len: " << samples.size() << std::endl;

    const int toSave = 100; // размер сохраняемого блока
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
        std::cout << "sample: ";
        for ( int k = 0; k < DIM; ++k )
            std::cout << sample[k] << " ";
        std::cout << std::endl;
        ++samples_counter;

        MPI_Send( &sample[0], DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent point " << endl;

        MPI_Send( &sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
        //cout << "(master) sent number of trajectory" << endl;

        ++sent;
    }

    std::vector<double> correlation_total(CORRELATION_FUNCTION_LENGTH);
    std::vector<double> correlation_package(CORRELATION_FUNCTION_LENGTH);

    int filesCounter = 0; // номер файла
    std::string filenameCorrelationFunction = "../output/eqcorr";

    double Volume = 4.0 / 3.0 * M_PI * pow( Rint_max * constants::ALU, 3);

    while( true )
    {
        std::cout << "(master) cycle beginning" << std::endl;

        if ( is_finished )
        {
            for ( int i = 1; i < world_size; i++ )
                MPI_Send( &is_finished, 1, MPI_INT, i, tags::EXIT_TAG, MPI_COMM_WORLD );

            break;
        }

        int msg;
        MPI_Recv( &msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        std::cout << "(master) received msg" << std::endl;

        source = status.MPI_SOURCE;
        if ( status.MPI_TAG == tags::TRAJECTORY_CUT_TAG )
        {
            std::cout << "(master) Received cutting trajectory tag!" << std::endl;
            if ( sent <= (int) samples.size() )
            {
                sample = samples[samples_counter];
                std::cout << "sample: ";
                for ( int k = 0; k < DIM; ++k )
                    std::cout << sample[k] << " ";
                std::cout << std::endl;
                ++samples_counter;

                MPI_Send( &sample[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
                MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
                //cout << "(master) Sent new point" << endl;

                continue;
            }
        }

        std::fill( correlation_package.begin(), correlation_package.end(), 0.0 );
        MPI_Recv( &correlation_package[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
        std::cout << "(master) received package" << std::endl;

        for ( size_t i = 0; i < CORRELATION_FUNCTION_LENGTH; i++ )
            correlation_total[i] += correlation_package[i];

        ++received;

        if ( (received % toSave == 0)  && (received != 0) )
        {
            std::string filename = filenameCorrelationFunction + std::to_string(filesCounter) + ".txt";
            ++filesCounter;

            std::cout << "(debug) Saving correlation file to " << filename << std::endl;
            save_correlation_function( correlation_total, filename,
                                       parameters.sampling_time * constants::ATU, Volume / received );
            // сохранили блок, обнулили элементы вектора
            std::fill( correlation_total.begin(), correlation_total.end(), 0.0 );
        }

        if ( received == (int) samples.size() )
        {
            std::cout << "(final) Exiting..." << std::endl;
            is_finished = true;
        }

        if ( sent < (int) samples.size() )
        {
            sample = samples[samples_counter];
            std::cout << "sample: ";
            for ( int k = 0; k < DIM; ++k )
                std::cout << sample[k] << " ";
            std::cout << std::endl;
            ++samples_counter;

            std::cout << "(master) sending initial condition" << std::endl;
            MPI_Send( &sample[0], DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
            MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
            std::cout << "(master) initial condition sent" << std::endl;

            ++sent;
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

    std::vector<double> dipx, dipy, dipz;
    std::vector<double> dipx_forward, dipy_forward, dipz_forward;
    std::vector<double> dipx_backward, dipy_backward, dipz_backward;
    std::vector<double> correlation_function;

    while ( true )
    {
        std::cout << "(slave) in the beginning of the cycle" << std::endl;

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
        std::cout << "(slave) received initial condition" << std::endl;

        if ( exit_status )
        {
            std::cout << "(" << world_rank << ") exit message is received!" << std::endl;
            break;
        }

        // начинаем расчет траектории
        trajectory.run_trajectory( syst );
        std::cout << "(slave) traversed trajectory" << std::endl;

        // если мы прошли предыдущий блок кода, значит траектория имеет допустимую длину.
        // копируем компоненты дипольного момента вдоль траектории в виде структуры данных vector<double>
        dipx_forward = trajectory.get_dipx();
        dipy_forward = trajectory.get_dipy();
        dipz_forward = trajectory.get_dipz();

        // после копирования освобождаем эти вектора внутри объекта trajectory
        trajectory.dump_dipoles( );

        // Прогоняем обратную траекторию
        trajectory.reverse_initial_conditions();

        trajectory.run_trajectory( syst );
        // собираем статус траектории. метод отправляет статус траектории мастеру
        // если траектория обрублена, то обнуляем буферы хранения диполя с траектории
        cut_trajectory = trajectory.report_trajectory_status( );
        if ( cut_trajectory )
        {
            trajectory.dump_dipoles();
            continue;
        }

        dipx_backward = trajectory.get_dipx();
        dipy_backward = trajectory.get_dipy();
        dipz_backward = trajectory.get_dipz();
        trajectory.dump_dipoles();

        // сливаем диполи с прямой и обратной траектории в одну
        merge(dipx, dipx_backward, dipx_forward);
        merge(dipy, dipy_backward, dipy_forward);
        merge(dipz, dipz_backward, dipz_forward);

        dipx_forward.clear();
        dipy_forward.clear();
        dipz_forward.clear();
        dipx_backward.clear();
        dipy_backward.clear();
        dipz_backward.clear();

        std::cout << "(" << world_rank << ") Processing " << trajectory.get_trajectory_counter()
                  << " trajectory. npoints = " << dipz_forward.size() << "; time = "
                  << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;

        calculate_physical_correlation(correlation_function, dipx, dipy, dipz);
        std::cout << "(slave) calculated physical correlation" << std::endl;

        dipx.clear();
        dipy.clear();
        dipz.clear();

        // Отправляем собранный массив корреляций мастер-процессу
        MPI_Send( &correlation_function[0], CORRELATION_FUNCTION_LENGTH, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
        correlation_function.clear();
        std::cout << "(slave) sent physical correlation" << std::endl;
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
