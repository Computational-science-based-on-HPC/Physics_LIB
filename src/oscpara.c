//
// Created by jghal on 6/16/2023.
//
/**
 * @file oscpara.c
 * @brief This file contains the implementation of the parallel versions of the oscillation simulation in 1D and 2D.
 *
 */
#include "../include/oscpara.h"
#include "../include/utils.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <fcntl.h>
#include <unistd.h>
#define NUM_THREADS 3
/**
 *  This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return
 */
int _simulate_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                         double Vo, double FI,
                                         double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("Invalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        max_amplitude = length;
        puts("Max Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }

    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    char *dir[2076];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    double buff[1000];
    if (world_rank < 1)
    {
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        coefficient_calc = -damping_coefficent / (2 * mass);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        int _sent_it_number = 0;

        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size))
                    _it_number_proc = (_it_number_all / (world_size));
                _sent_it_number += _it_number_proc;

                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        _it_number = _it_number_all - _sent_it_number;
    }

    if (world_rank > 0 && world_size > 1)
    {
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dir, sizeof(dir), MPI_CHAR, 0, MPI_COMM_WORLD);

    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(3);
    int count = 0;
    FILE *p_dis;
    char _file_name[2076];

    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    sprintf(_file_name, "%d_damped_os_parallel_v1_displacement_%d-%02d-%02d %02d:%02d:%02d.txt", world_rank, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    p_dis = fopen(_file_name, "w");

    for (int it = 0; it < _it_number; ++it)
    {

        if (RESULTS[0] == 0.000000 && RESULTS[1] == 0.000000 && RESULTS[2] == 0.000000)
        {
            _is_zero++;
            if (_is_zero > 2)
            {
                return 0;
            }
        }
        else
            _is_zero = 0;
        if (isnan(RESULTS[0]) || isnan(RESULTS[1]) || isnan(RESULTS[2]))
        {
            puts("Simulation Got a NaN Value.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a NaN Value Occurred");
            break;
        }
        else if (isinf(RESULTS[0]) || isinf(RESULTS[1]) || isinf(RESULTS[2]))
        {
            puts("Simulation Got a INF.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a INF Value Occurred");
            break;
        }
        if (count == 9999)
        {
            printf("count:%d rank:%d\n", count, world_rank);
            for (int i = 0; i < count; i++)
            {
                fprintf(p_dis, "%.6f\n", buff[i]);
            }
            count = 0;
        }
        t = step_size * (((double)it) + (world_rank * _it_number));

#pragma omp parallel sections
        {
#pragma omp section
            CALCULATIONS[0] = cos(W * t + FI);

#pragma omp section
            CALCULATIONS[1] = sin(W * t + FI);

#pragma omp section
            CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        }
#pragma omp parallel sections
        {
#pragma omp section
            RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];

#pragma omp section
            RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);

#pragma omp section
            RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
        }
        buff[count++] = RESULTS[0];
    }
    if (count > 0)
    {

        for (int i = 0; i < count; i++)
        {
            fprintf(p_dis, "%.6f\n", buff[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(p_dis);
    return 0;
}
/**
 *  This function calculate the execution time of simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @return
 */
int _execution_time_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                               double Vo, double FI,
                                               double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int initialized;
    MPI_Status status;

    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    time_t tim = time(NULL);
    double start_time, end_time;
    struct tm tm = *localtime(&tim);
    // char file_name[255];
    // sprintf(file_name, "DOP1%d-%02d-%02d %02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    int ret = 0;
    if (world_rank == 0)
    {
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("Started Simulation of Damped Oscillation Implementation Using MPI and OpenMP at %d-%02d-%02d %02d:%02d:%02d with Parametes:\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        printf("Amplitude: %f \nSpring Length: %f \nMass: %f \nGravity: %f \nStifeness: %f \nInitial_Acceleration: %f \nInitial_Velocity: %f \nFI_Const: %f \nTime_Limit: %f \nStep_Size(dt): %f \nDamping_coefficient: %f\nNumber of Processes: %d\n", max_amplitude, length, mass, gravity, k, Ao, Vo, FI, time_limit, step_size, damping_coefficent, world_size);
        puts("================================================================================");
        int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent, number_of_files, 0);
        if (validation == 0)
        {
            puts("\n\nInvalid Arguments is Given");
            ret = -1;
        }
        if (validation == -1)
        {
            max_amplitude = length;
            puts("\n\nMax Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
        }
        printf("Memory ===========================================================================\n");
        printmem();
        printf("\n================================================================================\nCPUs ===========================================================================\n\n");
        cpu_inf();
        printf("\n=================================================================================\n\n");

        tm = *localtime(&tim);
        start_time = MPI_Wtime();
        printf("Started Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank);
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        coefficient_calc = -damping_coefficent / (2 * mass);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        int _sent_it_number = 0;

        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                    _it_number_proc = (_it_number_all / (world_size - 1));
                _sent_it_number += _it_number_proc;
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }

        _it_number = _it_number_all - _sent_it_number;
        end_time = MPI_Wtime();
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("\nEnded Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank, end_time - start_time);
    }

    if (world_rank > 0 && world_size > 1)
    {
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ret, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (ret != 0)
    {
        return ret;
    }
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(3);
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nStarted Job 2 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank);
    start_time = MPI_Wtime();
    for (int it = 0; it < _it_number; ++it)
    {
        t = step_size * (((double)it) + (world_rank * _it_number));
#pragma omp parallel sections
        {
#pragma omp section
            CALCULATIONS[0] = cos(W * t + FI);

#pragma omp section
            CALCULATIONS[1] = sin(W * t + FI);

#pragma omp section
            CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        }

#pragma omp parallel sections
        {
#pragma omp section
            RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];

#pragma omp section
            RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);

#pragma omp section
            RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
        }
    }
    end_time = MPI_Wtime();
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nEnded Job 2 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank, end_time - start_time);

    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0)
    {
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("\n\nEnded Simulation of Damped Oscillation Implementation Using MPI and OpenMP at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    }
    return 0;
}
/**
 *  This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @return
 */
int _simulate_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent, int number_of_files)
{

    int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("Invalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        max_amplitude = length;
        puts("Max Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }

    int initialized;
    MPI_Initialized(&initialized);

    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    MPI_Status status;

    int world_size;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    double buff[1000];

    if (world_rank == 0)
    {
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        coefficient_calc = -damping_coefficent / (2 * mass);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        int _sent_it_number = 0;

        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                {
                    _it_number_proc = (_it_number_all / (world_size - 1));
                }
                _sent_it_number += _it_number_proc;
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        _it_number = _it_number_all - _sent_it_number;
    }

    MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (world_rank > 0 && world_size > 1)
    {
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    int count = 0;
    FILE *p_dis;
    char _file_name[2076];
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    sprintf(_file_name, "%damped_os_parallel_v2_displacement%d-%02d-%02d %02d:%02d:%02d.txt", world_rank, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    p_dis = fopen(_file_name, "w");

    for (int it = 0; it < _it_number; ++it)
    {

        if (RESULTS[0] == 0.000000 && RESULTS[1] == 0.000000 && RESULTS[2] == 0.000000)
        {
            _is_zero++;
            if (_is_zero > 2)
            {
                break;
            }
        }
        else
            _is_zero = 0;
        if (isnan(RESULTS[0]) || isnan(RESULTS[1]) || isnan(RESULTS[2]))
        {
            puts("Simulation Got a NaN Value.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a NaN Value Occurred");
            break;
        }
        else if (isinf(RESULTS[0]) || isinf(RESULTS[1]) || isinf(RESULTS[2]))
        {
            puts("Simulation Got a INF.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a INF Value Occurred");
            break;
        }
        if (count == 9999)
        {
            printf("count:%d rank:%d\n", count, world_rank);
            for (int i = 0; i < count; i++)
            {
                fprintf(p_dis, "%.6f\n", buff[i]);
            }
            count = 0;
        }
        t = step_size * (((double)it) + (world_rank * _it_number));
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
        buff[count++] = RESULTS[0];
    }
    if (count > 0)
    {

        for (int i = 0; i < count; i++)
        {
            fprintf(p_dis, "%.6f\n", buff[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(p_dis);
    return 0;
}
/**
 *  This function calculate execution time simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @return
 */
int _execution_time_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                           double Vo, double FI,
                                           double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }

    int world_size;
    int world_rank;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    int ret = 0;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    time_t tim = time(NULL);
    double start_time, end_time;
    struct tm tm = *localtime(&tim);
    if (world_rank == 0)
    {
        printf("Started Simulation of Damped Oscillation Implementation Using MPI at %d-%02d-%02d %02d:%02d:%02d with Parametes:\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        printf("Amplitude: %f \nSpring Length: %f \nMass: %f \nGravity: %f \nStifeness: %f \nInitial Acceleration: %f \nInitial Velocity: %f \nFI Const: %f \nTime Limit: %f \nStep_Size(dt): %f \nDamping coefficient: %f\nNumper of Processes", max_amplitude, length, mass, gravity, k, Ao, Vo, FI, time_limit, step_size, damping_coefficent, world_size);
        puts("================================================================================");
        printf("Memory ===========================================================================\n");
        printmem();
        printf("\n================================================================================\nCPUs ===========================================================================\n\n");
        cpu_inf();
        printf("\n=================================================================================\n\n");
        int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                    number_of_files, 0);
        if (validation == 0)
        {
            puts("\n\nInvalid Arguments is Given");
            ret = -1;
        }
        if (validation == -1)
        {
            max_amplitude = length;
            puts("\n\nMax Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
        }
        tim = time(NULL);
        tm = *localtime(&tim);
        start_time = MPI_Wtime();
        printf("Started Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank);

        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        coefficient_calc = -damping_coefficent / (2 * mass);
        number_of_files = _min_int(number_of_files, 3);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;

        _it_number_all = _round(time_limit / step_size);
        int _sent_it_number = 0;
        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                {
                    _it_number_proc = (_it_number_all / (world_size - 1));
                }
                _sent_it_number += _it_number_proc;
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        _it_number = _it_number_all - _sent_it_number;
        end_time = MPI_Wtime();
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("\nEnded Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank, end_time - start_time);
    }
    if (world_rank > 0 && world_size > 1)
    {
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ret, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (ret != 0)
    {
        return ret;
    }
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nStarted Job 2 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank);
    start_time = MPI_Wtime();

    for (int it = 0; it < _it_number; ++it)
    {
        t = step_size * (((double)it) + (world_rank * _it_number));
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
    }

    end_time = MPI_Wtime();
    // send execution times to rank 0 to get maximum
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nEnded Job 2 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, world_rank, end_time - start_time);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0)
    {
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("\n\nEnded Simulation of Damped Oscillation Implementation Using MPI at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    }
    return end_time - start_time;
}
/**
 *  This function calculate execution time simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using Openmp, The for loop is divided upon multiple threads.
 *  @note tried using sections pragma ended up applying more overhead in the code so performance decreased
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @param num_of_threads number of threads needed to execute the code
 * @return
 */
double
_execution_time_damped_os_parallel_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                       double Ao,
                                       double Vo, double FI,
                                       double time_limit, double step_size, double damping_coefficent,
                                       int number_of_files, int number_of_threads)
{
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    printf("Started Simulation of Damped Oscillation Implementation Using OpenMP at %d-%02d-%02d %02d:%02d:%02d with Parametes:\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    printf("Amplitude: %f \nSpring Length: %f \nMass: %f \nGravity: %f \nStifeness: %f \nInitial Acceleration: %f \nInitial Velocity: %f \nFI Const: %f \nTime Limit: %f \nStep Size(dt): %f \nDamping coefficient: %f\nNumber of Threads: %d\n", max_amplitude, length, mass, gravity, k, Ao, Vo, FI, time_limit, step_size, damping_coefficent, number_of_threads);
    puts("================================================================================");
    printf("Memory ===========================================================================\n");
    printmem();
    printf("\n================================================================================\nCPUs ===========================================================================\n\n");
    cpu_inf();
    printf("\n=================================================================================\n\n");
    int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("\n\nInvalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        max_amplitude = length;
        puts("\n\nMax Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    double Wo = sqrt(k / mass);
    double W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
    double RESULTS[3];
    RESULTS[0] = max_amplitude;
    RESULTS[1] = Vo + gravity * 0;
    RESULTS[2] = Ao + gravity * exp((-damping_coefficent / (2 * mass)) * 0);
    double CALCULATIONS[3];
    int num_steps = (int)(time_limit / step_size) + 10;
    printf("\nStarted Calculation at: %d-%02d-%02d %02d:%02d:%02d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    double start_time = omp_get_wtime();

#pragma omp parallel for num_threads(number_of_threads) private(CALCULATIONS, RESULTS)
    for (int i = 0; i <= num_steps; i++)
    {
        double t = i * step_size;
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
    }
    double end_time = omp_get_wtime();
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nEnded Calculation at: %d-%02d-%02d %02d:%02d:%02d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, end_time - start_time);
    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\n\nEnded Simulation of Damped Oscillation Implementation Using OpenMP at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    return end_time - start_time;
}