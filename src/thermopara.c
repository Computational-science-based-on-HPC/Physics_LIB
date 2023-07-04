#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/thermopara.h"
#include "../include/thermoutils.h"
#include "../include/utils.h"
#include <fcntl.h>
#include <unistd.h>

#define ll long long
#define THREADS 16
// #define M_PI 3.14159265358979323846264338327

// #ifdef __linux__
#include "mpi.h"
#include "omp.h"
// #endif

long double _get_value_1D_mpi(double time_step, double space_step, ll x, unsigned ll t, long long precision)
{
    long double sum = 0.0, exponential, spaceXTerm, coeff;
    long double x_real = x * space_step;
    long double t_real = t * time_step;

    for (ll k = 0; k < precision; k++)
    {
        exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
        spaceXTerm = sin((long double)(2 * k + 1) * M_PI * x_real / 2);
        coeff = 1 / (2 * k + 1);
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

long double _get_value_1D_openmp(double time_step, double space_step, ll x, unsigned ll t, long long precision)
{
    long double sum = 0.0, exponential, spaceXTerm, coeff;
    long double x_real = x * space_step;
    long double t_real = t * time_step;

#pragma omp parallel for num_threads(THREADS) schedule(static) shared(sum, x_real, t_real, precision) private(exponential, spaceXTerm, coeff)
    for (ll k = 0; k < precision; k++)
    {
        exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
        spaceXTerm = sin((long double)(2 * k + 1) * M_PI * x_real / 2);
        coeff = 1 / (2 * k + 1);
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

long double _get_value_2D_mpi(double time_step,
                         double length, double space_step_x, double width, double space_step_y,
                         ll x, ll y, unsigned ll t,
                         long long precision)
{
    long double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
    long double x_real = x * space_step_x;
    long double y_real = y * space_step_y;
    long double t_real = t * time_step;

    for (ll m = 1; m < precision; ++m)
    {
        for (ll n = 1; n < precision; ++n)
        {
            exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
            spaceXTerm = sin((long double)m * M_PI * x_real / length);
            spaceYTerm = sin((long double)n * M_PI * y_real / width);
            // Find Amn constant and multiply it with the sum
            coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
            sum += coeff * exponential * spaceXTerm * spaceYTerm;
        }
    }

    sum *= 200 / (M_PI * M_PI);
    return sum;
}

long double _get_value_2D_openmp(double time_step,
                            double length, double space_step_x, double width, double space_step_y,
                            ll x, ll y, unsigned ll t,
                            long long precision)
{
    long double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
    long double x_real = x * space_step_x;
    long double y_real = y * space_step_y;
    long double t_real = t * time_step;

#pragma omp parallel for num_threads(THREADS) schedule(static) shared(sum, x_real, y_real, t_real, precision) private(exponential, spaceXTerm, spaceYTerm, coeff)
    for (ll m = 1; m < precision; ++m)
    {
#pragma omp parallel for schedule(static)
        for (ll n = 1; n < precision; ++n)
        {
            exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
            spaceXTerm = sin((long double)m * M_PI * x_real / length);
            spaceYTerm = sin((long double)n * M_PI * y_real / width);
            // Find Amn constant and multiply it with the sum
            coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
            sum += coeff * exponential * spaceXTerm * spaceYTerm;
        }
    }
    sum *= 200 / (M_PI * M_PI);
    return sum;
}

int _simulate_heat_transfer_1D_MPI(double time_step, double time_limit, double space_step, long long precision)
{
    // MPI_Init(NULL, NULL);
    int initialized;
    double length = 10.0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;

    int my_rank;     // rank of process
    int processesNo; // number of process
    unsigned ll numTimePoint, numTimePointPerProcess, numTimePointRemProcess;
    ll numSpacePoint;
    double start_time, end_time;

    int checkRem;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        start_time = MPI_Wtime();
        numTimePoint = _cal_num_time(time_step, time_limit);
        numSpacePoint = _cal_num_space(length, space_step);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for remaining process
    }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePoint, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        unsigned ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        unsigned ll startIndex, i;

        MPI_Recv(&checkRem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        ll size;
        if (checkRem > 0)
        {
            size = numTimePointPerProcess + 1;
        }
        else
        {
            size = numTimePointPerProcess;
        }

        MPI_Recv(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        i = startIndex;
        unsigned ll endIndex = startIndex + size;

        char _file_name[2076];
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);

        sprintf(_file_name, "simulate_heat_transfer_1D_MPI_%d_%d-%02d-%02d_%02d-%02d-%02d.txt", my_rank, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        fptr1 = fopen(_file_name, "w");

        for (; i < endIndex; ++i)
        {
            for (ll x = 0; x <= numSpacePoint; ++x)
            {
                fprintf(fptr1, "%Lf ", _get_value_1D_mpi(time_step, space_step, x, i, precision));
            }
            fprintf(fptr1, "\n");
        }
    }
    fclose(fptr1);

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        printf("The time taken in MPI_1d With I/O is: %f\n", end_time - start_time);
    }
    return 0;
}

int _simulate_heat_transfer_1D_OPENMP(double time_step, double time_limit, double space_step, long long precision)
{
    clock_t start_time = clock();
    double length = 10.0;
    FILE *fptr;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_1D_OPENMP_%d-%02d-%02d_%02d-%02d-%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr = fopen(_file_name, "w");

    unsigned ll numTimePoint = _cal_num_time(time_step, time_limit);
    ll numSpacePoint = _cal_num_space(length, space_step);

    for (unsigned ll t = 0; t < numTimePoint; t++)
    {
        for (ll x = 0; x <= numSpacePoint; x++)
        {
            fprintf(fptr, "%Lf ", _get_value_1D_openmp(time_step, space_step, x, t, precision));
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_OPENMP_withFiles is: %f\n", execution_time);
    return 0;
}

int _simulate_heat_transfer_2D_MPI(double time_step, double time_limit,
                                   double space_step_x,
                                   double space_step_y,
                                   long long precision)
{
    int initialized;
    double length = 2.0, width = 2.0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;
    FILE *logFile;
    int my_rank;     // rank of process
    int processesNo; // number of process

    int checkRem;
    double start_time, end_time;

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    MPI_Barrier(MPI_COMM_WORLD);

    unsigned ll numTimePoint, numTimePointPerProcess, numTimePointRemProcess;
    ll numSpacePointX, numSpacePointY;


    if (my_rank == 0)
    {
        start_time = MPI_Wtime();
        numTimePoint = _cal_num_time(time_step, time_limit);
        numSpacePointX = _cal_num_space(length, space_step_x);
        numSpacePointY = _cal_num_space(width, space_step_y);


        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for remaining process
 }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointX, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointY, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        unsigned ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
        end_time = MPI_Wtime();
    }
    if (my_rank > 0)
    {
        unsigned ll startIndex, i;

        MPI_Recv(&checkRem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        ll size;
        if (checkRem > 0)
        {
            size = numTimePointPerProcess + 1;
        }
        else
        {
            size = numTimePointPerProcess;
        }

        MPI_Recv(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        i = startIndex;
        unsigned ll endIndex = startIndex + size;

        char _file_name[2076];
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        sprintf(_file_name, "simulate_heat_transfer_2D_MPI_%d_%d-%02d-%02d_%02d-%02d-%02d.txt", my_rank, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        fptr1 = fopen(_file_name, "w");

        for (; i < endIndex; ++i)
        {
            for (ll y = 0; y <= numSpacePointY; ++y)
            {
                for (ll x = 0; x <= numSpacePointX; ++x)
                {
                    fprintf(fptr1, "%Lf ", _get_value_2D_mpi(time_step, length, space_step_x, width, space_step_y, x, y, i, precision));
                }
                fprintf(fptr1, "\n");
            }
            fprintf(fptr1, "\n\n");
        }

        fclose(fptr1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        printf("The value of execution_time 2D_MPI_withFiles is: %f\n", execution_time);
    }
    
    return 0;
}

int _simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                      double space_step_x,
                                      double space_step_y,
                                      long long precision)
{

    clock_t start_time = clock();
    double length = 2.0, width = 2.0;
    FILE *fptr;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_2D_OPENMP_%d-%02d-%02d_%02d-%02d-%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr = fopen(_file_name, "w");

    unsigned ll numTimePoint;
    ll numSpacePointX, numSpacePointY;

    numTimePoint = _cal_num_time(time_step, time_limit);
    numSpacePointX = _cal_num_space(length, space_step_x);
    numSpacePointY = _cal_num_space(width, space_step_y);

    for (unsigned ll t = 0; t < numTimePoint; ++t)
    {
        for (ll y = 0; y <= numSpacePointY; ++y)
        {
            for (ll x = 0; x <= numSpacePointX; ++x)
            {
                fprintf(fptr, "%Lf ", _get_value_2D_openmp(time_step, length, space_step_x, width, space_step_y, x, y, t, precision));
            }
            fprintf(fptr, "\n");
        }
        fprintf(fptr, "\n\n");
    }

    fclose(fptr);
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_OPENMP_withFiles is: %f\n", execution_time);
    return 0;
}

double _execution_time_heat_transfer_1D_MPI(double time_step, double time_limit,
                                            double space_step,
                                            long long precision)
{
    FILE *logFile;
    int initialized;
    double length = 10.0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }

    char start_time_string[30], end_time_string[30];

    int my_rank;     // rank of process
    int processesNo; // number of process
    unsigned ll numTimePointPerProcess, numTimePointRemProcess, numTimePoint;
    ll  numSpacePoint;
    double start_time, end_time;

    int checkRem;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);

    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    if (my_rank == 0)
    {
        char _log_file_name[255];
        sprintf(_log_file_name, "Logs/Thermo Simulation MPI 1D/Thermo1D_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        logFile = fopen(_log_file_name, "w"); // Open the file in write mode
        if (logFile != NULL)
        {
            freopen(_log_file_name, "w", stdout);
        }
        tim = time(NULL);
        tm = *localtime(&tim);

        printf("Started Simulation of heat Equation 1D using MPI at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        printf("Number of processes: %d\n", processesNo);
        printf("Time step: %f\n", time_step);
        printf("Time limit: %f\n", time_limit);
        printf("Space step x: %f\n", space_step);
        printf("Precision: %lld\n", precision);
        printf("Length: %f\n", length);

        tim = time(NULL);
        tm = *localtime(&tim);



        start_time = MPI_Wtime();
        numTimePoint = _cal_num_time(time_step, time_limit);
        numSpacePoint = _cal_num_space(length, space_step);

        printf("Number of time points: %llu\n", numTimePoint);
        printf("Number of space points x: %lld\n", numSpacePoint);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process

        printf("Number of time points per process: %llu\n", numTimePointPerProcess);
        printf("Number of time points for remaining process: %llu\n", numTimePointRemProcess);
        puts("Memory ===========================================================================\n");
        printmemstream();
        puts("\n================================================================================\nCPUs ===========================================================================\n\n");
        cpu_inf_stream();
        puts("\n=================================================================================\n\n");

        printf("Started Job 1 at: %02d-%02d-%02d %02d:%02d:%02d in Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, my_rank);

    }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePoint, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        unsigned ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
        printf("\nEnded Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, my_rank, end_time - start_time);
    }
    else
    {
        unsigned ll startIndex, i;
        ll size;

        MPI_Recv(&checkRem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        if (checkRem > 0)
        {
            size = numTimePointPerProcess + 1;
        }
        else
        {
            size = numTimePointPerProcess;
        }

        MPI_Recv(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

//        printf("The value of startIndex is: %lld\n", startIndex);
        tim = time(NULL);
        tm = *localtime(&tim);
        start_time = MPI_Wtime();

        sprintf(start_time_string, "%04d-%02d-%02d %02d:%02d:%02d",
                tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                tm.tm_hour, tm.tm_min, tm.tm_sec);

        MPI_Send(start_time_string, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        i = startIndex;
        unsigned ll endIndex = startIndex + size;

        for (; i < endIndex; ++i)
        {
            for (ll x = 0; x <= numSpacePoint; ++x)
            {
                _get_value_1D_mpi(time_step, space_step, x, i, precision);
            }
        }

        tim = time(NULL);
        tm = *localtime(&tim);
        end_time = MPI_Wtime();

        sprintf(end_time_string, "%04d-%02d-%02d %02d:%02d:%02d",
                tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                tm.tm_hour, tm.tm_min, tm.tm_sec);

        MPI_Send(end_time_string, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        double execution = end_time - start_time ;
        MPI_Send(&execution, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double execution_time;
    if (my_rank == 0)
    {
        double execution_recv;
        for (int i = 1; i < processesNo; i++) {
            char rank_start_time[30], rank_end_time[30];
            MPI_Recv(rank_start_time, 30, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(rank_end_time, 30, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&execution_recv, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);

            printf("Rank %d start time: %s, end time: %s Execution Time: %f sec.\n", i, rank_start_time, rank_end_time, execution_recv);
        }
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;

        tim = time(NULL);
        tm = *localtime(&tim);
        printf("Finished the execution of heat Equation 2D using MPI at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        fclose(logFile);
    }

    return execution_time;
}

double _execution_time_heat_transfer_1D_OPENMP(double time_step, double time_limit,
                                               double space_step,
                                               long long precision)
{

    FILE *file;
    double length = 10.0;
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    char _log_file_name[255];
    sprintf(_log_file_name, "Logs/Thermo Simulation execution openmp 1D/Thermo1D_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    file = fopen(_log_file_name, "w"); // Open the file in write mode
    if (file != NULL)
    {
        freopen(_log_file_name, "w", stdout);
    }

    printf("Started Simulation of heat Equation 1D using OPENMP at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    printf("Time step: %f\n", time_step);
    printf("Time limit: %f\n", time_limit);
    printf("Space step x: %f\n", space_step);
    printf("Precision: %lld\n", precision);
    printf("Length: %f\n", length);

    unsigned ll numTimePoint = _cal_num_time(time_step, time_limit);
    ll numSpacePoint = _cal_num_space(length, space_step);

    printf("Number of time points: %llu\n", numTimePoint);
    printf("Number of space points x: %lld\n", numSpacePoint);
    printf("Memory ===========================================================================\n");
    printmemstream();
    printf("\n================================================================================\nCPUs ===========================================================================\n\n");
    cpu_inf_stream();
    printf("\n=================================================================================\n\n");

    printf("\nStarted Calculation at: %d-%02d-%02d_%02d-%02d-%02d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    double start_time2 = omp_get_wtime();
    for (unsigned ll t = 0; t < numTimePoint; t++)
    {
        for (ll x = 0; x <= numSpacePoint; x++)
        {
            _get_value_1D_openmp(time_step, space_step, x, t, precision);
        }
    }
    double end_time2 = omp_get_wtime();

    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nEnded Calculation at: %d-%02d-%02d_%02d-%02d-%02d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, end_time2 - start_time2);
    tim = time(NULL);
    tm = *localtime(&tim);
    puts("\n================================================================================\n");
    printf("Ended Simulation of Thermo energy 1D at %d-%02d-%02d_%02d-%02d-%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

    return end_time2 - start_time2; // execution_time;
}

double _execution_time_heat_transfer_2D_MPI(double time_step, double time_limit,
                                            double space_step_x,
                                            double space_step_y,
                                            long long precision)
{

    FILE *logFile;
    int initialized;
    double length = 2.0, width = 2.0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    // MPI_Init(NULL, NULL);
    int my_rank;     // rank of process
    int processesNo; // number of process

    int checkRem;
    double start_time, end_time;
    char start_time_string[30], end_time_string[30];
    
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);

    unsigned ll numTimePoint, numTimePointPerProcess, numTimePointRemProcess;
    ll numSpacePointX, numSpacePointY;

    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    if (my_rank == 0)
    {
        char _log_file_name[255];
        sprintf(_log_file_name, "Logs/Thermo Simulation mpi 2D/Thermo2D_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        logFile = fopen(_log_file_name, "w"); // Open the file in write mode
        if (logFile != NULL)
        {
            freopen(_log_file_name, "w", stdout);
        }
        tim = time(NULL);
        tm = *localtime(&tim);

        printf("Started Simulation of heat Equation 2D using MPI at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        printf("Number of processes: %d\n", processesNo);
        printf("Time step: %f\n", time_step);
        printf("Time limit: %f\n", time_limit);
        printf("Space step x: %f\n", space_step_x);
        printf("Space step y: %f\n", space_step_y);
        printf("Precision: %lld\n", precision);
        printf("Length: %f\n", length);
        printf("Width: %f\n", width);

        start_time = MPI_Wtime();
        tim = time(NULL);
        tm = *localtime(&tim);

        numTimePoint = _cal_num_time(time_step, time_limit);
        numSpacePointX = _cal_num_space(length, space_step_x);
        numSpacePointY = _cal_num_space(width, space_step_y);

        printf("Number of time points: %llu\n", numTimePoint);
        printf("Number of space points x: %lld\n", numSpacePointX);
        printf("Number of space points y: %lld\n", numSpacePointY);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for remaining process

        printf("Number of time points per process: %llu\n", numTimePointPerProcess);
        printf("Number of time points for remaining process: %llu\n", numTimePointRemProcess);
        puts("Memory ===========================================================================\n");
        printmemstream();
        puts("\n================================================================================\nCPUs ===========================================================================\n\n");
        cpu_inf_stream();
        puts("\n=================================================================================\n\n");

        printf("Started Job 1 at: %02d-%02d-%02d %02d:%02d:%02d in Processor: %s Rank: %d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, my_rank);
    }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointX, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointY, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        unsigned ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
        end_time = MPI_Wtime();
        tim = time(NULL);
        tm = *localtime(&tim);
        printf("\nEnded Job 1 at: %d-%02d-%02d %02d:%02d:%02d Processor: %s Rank: %d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, processor_name, my_rank, end_time - start_time);
    }
    if (my_rank > 0)
    {
        unsigned ll startIndex, i;
        start_time = MPI_Wtime();

        MPI_Recv(&checkRem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        ll size;
        if (checkRem > 0)
        {
            size = numTimePointPerProcess + 1;
        }
        else
        {
            size = numTimePointPerProcess;
        }

        MPI_Recv(&startIndex, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        tim = time(NULL);
        tm = *localtime(&tim);
        start_time = MPI_Wtime();

        sprintf(start_time_string, "%04d-%02d-%02d %02d:%02d:%02d",
        tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
        tm.tm_hour, tm.tm_min, tm.tm_sec);

        MPI_Send(start_time_string, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        i = startIndex;
        unsigned ll endIndex = startIndex + size;

        for (; i < endIndex; ++i)
        {
            for (ll y = 0; y <= numSpacePointY; ++y)
            {
                for (ll x = 0; x <= numSpacePointX; ++x)
                {
                    _get_value_2D_mpi(time_step, length, space_step_x, width, space_step_y, x, y, i, precision);
                }
            }
        }

        tim = time(NULL);
        tm = *localtime(&tim);
        end_time = MPI_Wtime();

        sprintf(end_time_string, "%04d-%02d-%02d %02d:%02d:%02d",
        tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
        tm.tm_hour, tm.tm_min, tm.tm_sec);

        MPI_Send(end_time_string, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        double execution = end_time - start_time ;
        MPI_Send(&execution, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double execution_time;
    if (my_rank == 0)
    {
        double execution_recv;
        for (int i = 1; i < processesNo; i++) {
            char rank_start_time[30], rank_end_time[30];
            MPI_Recv(rank_start_time, 30, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(rank_end_time, 30, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&execution_recv, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);

            printf("Rank %d start time: %s, end time: %s Execution Time: %f sec.\n", i, rank_start_time, rank_end_time, execution_recv);
        }
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;

        tim = time(NULL);
        tm = *localtime(&tim);
        printf("Finished the execution of heat Equation 2D using MPI at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        fclose(logFile);
    }
    return execution_time;
}

double
_execution_time_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                        double space_step_x,
                                        double space_step_y,
                                        long long precision)
{

    //    clock_t start_time = clock();
    double length = 2.0, width = 2.0;

    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    char _log_file_name[255];
    FILE *file;
    sprintf(_log_file_name, "Logs/Thermo Simulation execution openmp 2D/Thermo2D_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    file = fopen(_log_file_name, "w"); // Open the file in write mode
    if (file != NULL)
    {
        freopen(_log_file_name, "w", stdout);
    }

    printf("Started Simulation of heat Equation 2D using OPENMP at %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    printf("Time step: %f\n", time_step);
    printf("Time limit: %f\n", time_limit);
    printf("Space step x: %f\n", space_step_x);
    printf("Space step y: %f\n", space_step_y);
    printf("Precision: %lld\n", precision);
    printf("Length: %f\n", length);
    printf("Width: %f\n", width);

    unsigned ll numTimePoint;
    ll numSpacePointX, numSpacePointY;

    numTimePoint = _cal_num_time(time_step, time_limit);
    numSpacePointX = _cal_num_space(length, space_step_x);
    numSpacePointY = _cal_num_space(width, space_step_y);
    printf("Number of time points: %llu\n", numTimePoint);
    printf("Number of space points x: %lld\n", numSpacePointX);
    printf("Number of space points y: %lld\n", numSpacePointY);
    printf("Memory ===========================================================================\n");
    printmemstream();
    printf("\n================================================================================\nCPUs ===========================================================================\n\n");
    cpu_inf_stream();
    printf("\n=================================================================================\n\n");

    printf("\nStarted Calculation at: %d-%02d-%02d %02d:%02d:%02d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    double start_time = omp_get_wtime();

    for (unsigned ll t = 0; t < numTimePoint; ++t)
    {
        for (ll y = 0; y <= numSpacePointY; ++y)
        {
            for (ll x = 0; x <= numSpacePointX; ++x)
            {
                _get_value_2D_openmp(time_step, length, space_step_x, width, space_step_y, x, y, t, precision);
            }
        }
    }

    double end_time = omp_get_wtime();
    //    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    //    printf("The value of execution_time 2D_OPENMP_without_Files is: %f\n", end_time - start_time);

    tim = time(NULL);
    tm = *localtime(&tim);
    printf("\nEnded Calculation at: %d-%02d-%02d %02d:%02d:%02d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, end_time - start_time);
    tim = time(NULL);
    tm = *localtime(&tim);
    puts("\n================================================================================\n");
    printf("Ended Simulation of Thermo energy 2D at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

    return end_time - start_time; // execution_time;
}
