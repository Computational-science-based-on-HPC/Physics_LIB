//
// Created by jghal on 6/16/2023.
//
#include "oscpara.h"
#include "utils.h"
#include "stdio.h"
#include "math.h"
#include "time.h"

#ifdef MPI_INCLUDE
#include "mpi.h"
#include "omp.h"
#endif
#define NUM_THREADS 3

const char *FILES[] = {"displacement.txt", "velocity.txt", "acceleration.txt"};

int
_simulate_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                         double Vo, double FI,
                                         double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    MPI_Init(NULL, NULL);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;
    MPI_Offset _d_offset;
    MPI_File fh;
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    if (world_rank == 0)
    {
        int validation = _valid_osc(max_amplitude, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
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
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        RESULTS[3];
        coefficient_calc = -damping_coefficent / (2 * mass);
        CALCULATIONS[3];
        number_of_files = _min_int(number_of_files, 3);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                    _it_number_proc = (_it_number_all / (world_size - 1));
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            MPI_Bcast(&number_of_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    if (world_rank > 0 && world_size > 1)
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    else
    {
        _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
        if (0 > _it_number_all % (world_size - 1))
            _it_number_proc = (_it_number_all / (world_size - 1));
        _it_number = _it_number_proc;
    }
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(NUM_THREADS);
    int count = 0;
    for (int it = world_rank * _it_number; it < _it_number_all; ++it)
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
        t = step_size * ((double)it);
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
    MPI_Finalize();
}
double
_execution_time_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                         double Vo, double FI,
                                         double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    MPI_Init(NULL, NULL);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;
    MPI_Offset _d_offset;
    MPI_File fh;
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    if (world_rank == 0)
    {
        int validation = _valid_osc(max_amplitude, 0,length, mass, gravity, k, time_limit, step_size, damping_coefficent,
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
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        RESULTS[3];
        coefficient_calc = -damping_coefficent / (2 * mass);
        CALCULATIONS[3];
        number_of_files = _min_int(number_of_files, 3);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                    _it_number_proc = (_it_number_all / (world_size - 1));
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            MPI_Bcast(&number_of_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    if (world_rank > 0 && world_size > 1)
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    else
    {
        _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
        if (0 > _it_number_all % (world_size - 1))
            _it_number_proc = (_it_number_all / (world_size - 1));
        _it_number = _it_number_proc;
    }
    clock_t start_time = clock();
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(NUM_THREADS);
    int count = 0;
    for (int it = world_rank * _it_number; it < _it_number_all; ++it)
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
        t = step_size * ((double)it);
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
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // TODO if not 0 send to rank 0 execution time
    MPI_Finalize();
}

int
_simulate_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    MPI_Init(NULL, NULL);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;
    MPI_Offset _d_offset;
    MPI_File fh;
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    if (world_rank == 0)
    {
        int validation = _valid_osc(max_amplitude, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
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
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        RESULTS[3];
        coefficient_calc = -damping_coefficent / (2 * mass);
        CALCULATIONS[3];
        number_of_files = _min_int(number_of_files, 3);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                    _it_number_proc = (_it_number_all / (world_size - 1));
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            MPI_Bcast(&number_of_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    if (world_rank > 0 && world_size > 1)
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    else
    {
        _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
        if (0 > _it_number_all % (world_size - 1))
            _it_number_proc = (_it_number_all / (world_size - 1));
        _it_number = _it_number_proc;
    }
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(NUM_THREADS);
    int count = 0;
    for (int it = world_rank * _it_number; it < _it_number_all; ++it)
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
        t = step_size * ((double)it);
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
    }
    MPI_Finalize();
}
double
_execution_time_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                       double Vo, double FI,
                                       double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    MPI_Init(NULL, NULL);
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status status;
    MPI_Offset _d_offset;
    MPI_File fh;
    double W;
    double RESULTS[3];
    double coefficient_calc;
    double CALCULATIONS[3];
    int _it_number_proc;
    int _it_number;
    int _it_number_all;
    double t;
    short int _is_zero = 0;
    if (world_rank == 0)
    {
        int validation = _valid_osc(max_amplitude, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
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
        double Wo = sqrt(k / mass);
        W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
        RESULTS[3];
        coefficient_calc = -damping_coefficent / (2 * mass);
        CALCULATIONS[3];
        number_of_files = _min_int(number_of_files, 3);
        RESULTS[0] = max_amplitude;
        RESULTS[1] = Vo;
        RESULTS[2] = Ao + gravity;
        _it_number_all = _round(time_limit / step_size);
        if (world_size > 1)
        {
            _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
            for (int i = 1; i < world_size; i++)
            {
                if (i > _it_number_all % (world_size - 1))
                    _it_number_proc = (_it_number_all / (world_size - 1));
                MPI_Send(&_it_number_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            MPI_Bcast(&number_of_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&W, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&coefficient_calc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    if (world_rank > 0 && world_size > 1)
        MPI_Recv(&_it_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    else
    {
        _it_number_proc = (_it_number_all / (world_size - 1)) + 1;
        if (0 > _it_number_all % (world_size - 1))
            _it_number_proc = (_it_number_all / (world_size - 1));
        _it_number = _it_number_proc;
    }
    clock_t start_time = clock();
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(NUM_THREADS);
    int count = 0;
    for (int it = world_rank * _it_number; it < _it_number_all; ++it)
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
        t = step_size * ((double)it);
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
    }
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // TODO if not zero send execution time
    MPI_Finalize();
}