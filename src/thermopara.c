#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/thermopara.h"
#include "../include/thermoutils.h"
#define ll long long
#define THREADS 8
// #define M_PI 3.14159265358979323846264338327

// #ifdef __linux__
#include "mpi.h"
#include "omp.h"
// #endif

double _get_value_1D_mpi(double time_step, double space_step, double x, double t, int precision)
{
    double sum = 0.0, exponential, spaceXTerm, coeff;
    double x_real = x * space_step;
    double t_real = t * time_step;

    for (int k = 0; k < precision; k++)
    {
        exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
        spaceXTerm = sin((double)(2 * k + 1) * M_PI * x_real / 2);
        coeff = 1 / (2 * k + 1);
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

double _get_value_1D_openmp_V1(double time_step, double space_step, double x, double t, int precision){
    double sum = 0.0, exponential, spaceXTerm, coeff;
    double x_real = x * space_step;
    double t_real = t * time_step;

    #pragma omp parallel for num_threads(THREADS) schedule(static) shared(sum, x_real, t_real, precision) private(exponential, spaceXTerm, coeff)
    for (int k = 0; k < precision; k++)
    {
        exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
        spaceXTerm = sin((double)(2 * k + 1) * M_PI * x_real / 2);
        coeff = 1 / (2 * k + 1);
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

double _get_value_1D_openmp_V2(double time_step, double space_step, double x, double t, int precision){
    double sum = 0.0, exponential, spaceXTerm, coeff;
    double x_real = x * space_step;
    double t_real = t * time_step;

    #pragma omp parallel for num_threads(THREADS) schedule(static) shared(sum, x_real, t_real, precision) private(exponential, spaceXTerm, coeff)
    for (int k = 0; k < precision; k++)
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
            #pragma omp section
            spaceXTerm = sin((double)(2 * k + 1) * M_PI * x_real / 2);
            #pragma omp section
            coeff = 1 / (2 * k + 1);
        }
        
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

 double _get_value_2D_mpi(double time_step,
                          double length, double space_step_x, double width, double space_step_y,
                          int x, int y, int t,
                          int precision){
     double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
     double x_real = x * space_step_x;
     double y_real = y * space_step_y;
     double t_real = t * time_step;

     for (ll m = 1; m < precision; ++m)
     {
         for (ll n = 1; n < precision; ++n)
         {
             exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
             spaceXTerm = sin((double)m * M_PI * x_real / length);
             spaceYTerm = sin((double)n * M_PI * y_real / width);
             // Find Amn constant and multiply it with the sum
             coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
             sum += coeff * exponential * spaceXTerm * spaceYTerm;
         }
     }

     sum *= 200 / (M_PI * M_PI);
     return sum;
 }

 double _get_value_2D_openmp(double time_step,
                          double length, double space_step_x, double width, double space_step_y,
                          int x, int y, int t,
                          int precision){
     double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
     double x_real = x * space_step_x;
     double y_real = y * space_step_y;
     double t_real = t * time_step;

     #pragma omp parallel for num_threads(THREADS) schedule(static) shared(sum, x_real, y_real, t_real, precision) private(exponential, spaceXTerm, spaceYTerm, coeff)
     for (ll m = 1; m < precision; ++m)
     {
         #pragma omp parallel for schedule(static)
         for (ll n = 1; n < precision; ++n)
         {
             exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
             spaceXTerm = sin((double)m * M_PI * x_real / length);
             spaceYTerm = sin((double)n * M_PI * y_real / width);
             // Find Amn constant and multiply it with the sum
             coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
             sum += coeff * exponential * spaceXTerm * spaceYTerm;
         }
     }
     sum *= 200 / (M_PI * M_PI);
     return sum;
 }

 double _get_value_2D_openmp_v2(double time_step,
                          double length, double space_step_x, double width, double space_step_y,
                          int x, int y, int t,
                          int precision){
     double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
    double x_real = x * space_step_x;
    double y_real = y * space_step_y;
    double t_real = t * time_step;

     #pragma omp parallel for num_threads(THREADS) shared(sum, x_real, y_real, t_real, precision) private(exponential, spaceXTerm, spaceYTerm, coeff)
     for (ll m = 1; m < precision; ++m)
     {
         // #pragma omp parallel for schedule(static)
         for (ll n = 1; n < precision; ++n)
         {
             #pragma omp parallel sections
         {
             #pragma omp section
             exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
             #pragma omp section
             spaceXTerm = sin((double)m * M_PI * x_real / length);
             #pragma omp section
             spaceYTerm = sin((double)n * M_PI * y_real / width);
             #pragma omp section
             // Find Amn constant and multiply it with the sum
             coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);


         }
             sum += coeff * exponential * spaceXTerm * spaceYTerm;
         }
     }
     sum *= 200 / (M_PI * M_PI);
     return sum;
 }


int _simulate_heat_transfer_1D_MPI(double time_step, double time_limit, double length, double space_step, int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
     MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    fptr1 = fopen("1_1D_MPI.txt", "w");
    fptr2 = fopen("2_1D_MPI.txt", "w");
    fptr3 = fopen("3_1D_MPI.txt", "w");
    int my_rank;     // rank of process
    int processesNo; // number of process
    ll numTimePointPerProcess, numTimePointRemProcess, numTimePoint, numSpacePoint;
    double start_time, end_time;

    int checkRem;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);


    if(my_rank == 0){
        start_time = MPI_Wtime();
        numTimePoint= _cal_num_time(time_step, time_limit);
        numSpacePoint= _cal_num_space(length, space_step);

        // printf("The value of numTimePoint is: %lld\n", numTimePoint);
        // printf("The value of numSpacePoint is: %lld\n", numSpacePoint);
        // printf("The value of time_step is: %f\n", time_step);
        // printf("The value of length is: %f\n", length);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process
    }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePoint, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        ll startIndex;
        ll i;

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

        MPI_Recv(&startIndex, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        printf("The value of startIndex is: %lld\n", startIndex);

        i = startIndex;
        ll endIndex = startIndex + size;

        for (; i < endIndex; ++i)
        {
            for (ll x = 0; x < numSpacePoint; ++x)
            {
                if (my_rank == 1)
                {
                    fprintf(fptr1, "%f ", _get_value_1D_mpi(time_step, space_step, x, i, precision));
                }
                else if (my_rank == 2)
                {
                    fprintf(fptr2, "%f ", _get_value_1D_mpi(time_step, space_step, x, i, precision));
                }
                else if (my_rank == 3)
                {
                    fprintf(fptr3, "%f ",_get_value_1D_mpi(time_step, space_step, x, i, precision));
                }
            }
            if (my_rank == 1)
            {
                fprintf(fptr1, "\n");
            }
            else if (my_rank == 2)
            {
                fprintf(fptr2, "\n");
            }
            else if (my_rank == 3)
            {
                fprintf(fptr3, "\n");
            }
        }
    }
    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        printf("The time taken in MPI_1d With I/O is: %f\n", end_time - start_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
    return 0;
}

int _simulate_heat_transfer_1D_OPENMP(double time_step, double time_limit, double length, double space_step, int precision){
    clock_t start_time=clock();
    FILE *fptr;
    fptr = fopen("1D_OPENMP_V1.txt", "w");

    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePoint= _cal_num_space(length, space_step);

    for (int t = 0; t < numTimePoint; t++) {
        for (int x = 0; x < numSpacePoint; x++) {
            fprintf(fptr, "%f ", _get_value_1D_openmp_V1(time_step, space_step, x, t, precision));
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_OPENMP_withFiles is: %f\n",execution_time);
    return 0;
}

int
_simulate_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit, double length, double space_step, int precision){
    clock_t start_time=clock();
    FILE *fptr;
    fptr = fopen("1D_OPENMP_V2.txt", "w");

    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePoint= _cal_num_space(length, space_step);

    for (int t = 0; t < numTimePoint; t++) {
        for (int x = 0; x < numSpacePoint; x++) {
            fprintf(fptr, "%f ", _get_value_1D_openmp_V2(time_step, space_step, x, t, precision));
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_OPENMP_V2_withFiles is: %f\n",execution_time);
    return 0;

}


 int
 _simulate_heat_transfer_2D_MPI(double time_step, double time_limit,
                                double length, double space_step_x,
                                double width, double space_step_y,
                                int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
     MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    fptr1 = fopen("1_2D_MPI.txt", "w");
    fptr2 = fopen("2_2D_MPI.txt", "w");
    fptr3 = fopen("3_2D_MPI.txt", "w");

    int my_rank;     // rank of process
    int processesNo; // number of process

    int checkRem;
    double start_time, end_time;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);
     

    ll numTimePoint;
    ll numSpacePointX;
    ll numSpacePointY;
    ll numTimePointPerProcess;
    ll numTimePointRemProcess;

    //_cal_num_time(time_param, &numTimePoint);
    if(my_rank == 0){
        start_time = MPI_Wtime();
        numTimePoint= _cal_num_time(time_step, time_limit);
        numSpacePointX= _cal_num_space(length, space_step_x);
        numSpacePointY= _cal_num_space(width, space_step_y);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process
    }


    MPI_Bcast(&numTimePointPerProcess, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointX, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointY, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        ll startIndex;
        ll i;

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

        MPI_Recv(&startIndex, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        printf("The value of startIndex is: %lld\n", startIndex);

        i = startIndex;
        ll endIndex = startIndex + size;


        for (; i < endIndex; ++i)
        {
            for (ll y = 0; y <= numSpacePointY; ++y)
            {
                for (ll x = 0; x <= numSpacePointX; ++x)
                {
                    if(my_rank == 1){
                        fprintf(fptr1, "%f ", _get_value_2D_mpi(time_step, length, space_step_x, width, space_step_y, x, y, i, precision));
                    }else if(my_rank == 2){
                        fprintf(fptr2, "%f ", _get_value_2D_mpi(time_step, length, space_step_x, width, space_step_y, x, y, i, precision));
                    }else if(my_rank == 3){
                        fprintf(fptr3, "%f ", _get_value_2D_mpi(time_step, length, space_step_x, width, space_step_y, x, y, i, precision));
                    }
                
                }
                if(my_rank == 1){
                    fprintf(fptr1, "\n");
                }else if(my_rank == 2){
                    fprintf(fptr2, "\n");
                }else if(my_rank == 3){
                    fprintf(fptr3, "\n");
                }
                // fprintf(fptr, "\n");
            }
            if(my_rank == 1){
                fprintf(fptr1, "\n\n");
            }else if(my_rank == 2){
                fprintf(fptr2, "\n\n");
            }else if(my_rank == 3){
                fprintf(fptr3, "\n\n");
            }
            // fprintf(fptr, "\n\n");
        }
    }

    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        printf("The time taken in MPI_2d With I/O is: %f\n", end_time - start_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
            
    return 0;
 }

 int
 _simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double length, double space_step_x,
                                   double width, double space_step_y,
                                   int precision){

    clock_t start_time=clock();
    FILE *fptr;
    fptr = fopen("2D_OPENMP.txt", "w");

    ll numTimePoint;
    ll numSpacePointX;
    ll numSpacePointY;

    numTimePoint= _cal_num_time(time_step, time_limit);
    numSpacePointX= _cal_num_space(length, space_step_x);
    numSpacePointY= _cal_num_space(width, space_step_y);

    for (ll t = 0; t < numTimePoint; ++t) {
        for (ll y = 0; y <= numSpacePointY; ++y) {
            for (ll x = 0; x <= numSpacePointX; ++x) {
                fprintf(fptr, "%f ", _get_value_2D_openmp(time_step, length, space_step_x, width, space_step_y, x, y, t, precision));
            }
            fprintf(fptr, "\n");
        }
        fprintf(fptr, "\n\n");
    }

    fclose(fptr);
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_OPENMP_withFiles is: %f\n",execution_time);
    return 0;

 }

 int
 _simulate_heat_transfer_2D_OPENMP_V2(double time_step, double time_limit,
                                      double length, double space_step_x,
                                      double width, double space_step_y,
                                      int precision){

    clock_t start_time=clock();
    FILE *fptr;
    fptr = fopen("2D_OPENMP_V2.txt", "w");



    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePointX= _cal_num_space(length, space_step_x);
    ll numSpacePointY= _cal_num_space(width, space_step_y);
    for (ll t = 0; t < numTimePoint; ++t) {
        for (ll y = 0; y <= numSpacePointY; ++y) {
            for (ll x = 0; x <= numSpacePointX; ++x) {
                fprintf(fptr, "%f ", _get_value_2D_openmp_v2(time_step, length, space_step_x, width, space_step_y, x, y, t, precision));
            }
            fprintf(fptr, "\n");
        }
        fprintf(fptr, "\n\n");
    }

    fclose(fptr);
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_OPENMP_V2_withFiles is: %f\n",execution_time);
    return 0;

 }

int _execution_time_heat_transfer_1D_MPI(double time_step, double time_limit,
                                     double length, double space_step,
                                     int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
     MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }

    int my_rank;     // rank of process
    int processesNo; // number of process
    ll numTimePointPerProcess, numTimePointRemProcess, numTimePoint, numSpacePoint;
    double start_time, end_time;

    int checkRem;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        start_time = MPI_Wtime();
        numTimePoint= _cal_num_time(time_step, time_limit);
        numSpacePoint= _cal_num_space(length, space_step);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process
    }

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePoint, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        ll startIndex;
        ll i;

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

        MPI_Recv(&startIndex, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        printf("The value of startIndex is: %lld\n", startIndex);

        i = startIndex;
        ll endIndex = startIndex + size;

        for (; i < endIndex; ++i)
        {
            for (ll x = 0; x < numSpacePoint; ++x)
            {
                _get_value_1D_mpi(time_step, space_step, x, i, precision);
            }

        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        printf("The time taken in MPI_1d Without I/O is: %f\n", end_time - start_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
        
    return 0;
}
int _execution_time_heat_transfer_1D_OPENMP(double time_step, double time_limit,
                                        double length, double space_step,
                                        int precision){

    clock_t start_time=clock();

    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePoint= _cal_num_space(length, space_step);

    for (int t = 0; t < numTimePoint; t++) {
        for (int x = 0; x < numSpacePoint; x++) {
            _get_value_1D_openmp_V1(time_step, space_step, x, t, precision);
        }
    }
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_OPENMP is: %f\n",execution_time);
    return 0;
}

int _execution_time_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit,
                                           double length, double space_step,
                                           int precision){

    clock_t start_time=clock();

    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePoint= _cal_num_space(length, space_step);

    for (int t = 0; t < numTimePoint; t++) {
        for (int x = 0; x < numSpacePoint; x++) {
            _get_value_1D_openmp_V2(time_step, space_step, x, t, precision);
        }

    }
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_OPENMP_V2 is: %f\n",execution_time);
    return 0;

}
int _execution_time_heat_transfer_2D_MPI(double time_step, double time_limit,
                                     double length, double space_step_x,
                                     double width, double space_step_y,
                                     int precision){
    
    int initialized, finalized;
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
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);
    MPI_Barrier(MPI_COMM_WORLD);

    ll numTimePoint;
    ll numSpacePointX;
    ll numSpacePointY;
    ll numTimePointPerProcess;
    ll numTimePointRemProcess;

    //_cal_num_time(time_param, &numTimePoint);
    if(my_rank == 0){
        start_time = MPI_Wtime();
        numTimePoint= _cal_num_time(time_step, time_limit);
        numSpacePointX= _cal_num_space(length, space_step_x);
        numSpacePointY= _cal_num_space(width, space_step_y);

        numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
        numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process
    }


    MPI_Bcast(&numTimePointPerProcess, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointX, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePointY, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_step_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_limit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        ll startIndex;
        ll i;

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

        MPI_Recv(&startIndex, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        printf("The value of startIndex is: %lld\n", startIndex);

        i = startIndex;
        ll endIndex = startIndex + size;


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
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        printf("The time taken in MPI_2d Without I/O is: %f\n", end_time - start_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
    
    return 0;

}

 int
 _execution_time_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double length, double space_step_x,
                                   double width, double space_step_y,
                                   int precision){

    clock_t start_time=clock();
    ll numTimePoint;
    ll numSpacePointX;
    ll numSpacePointY;

    numTimePoint= _cal_num_time(time_step, time_limit);
    numSpacePointX= _cal_num_space(length, space_step_x);
    numSpacePointY= _cal_num_space(width, space_step_y);

    for (ll t = 0; t < numTimePoint; ++t) {
        for (ll y = 0; y <= numSpacePointY; ++y) {
            for (ll x = 0; x <= numSpacePointX; ++x) {
                _get_value_2D_openmp(time_step, length, space_step_x, width, space_step_y, x, y, t, precision);
            }
        }
    }

    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_OPENMP_without_Files is: %f\n",execution_time);
    return 0;

 }

 int
 _execution_time_heat_transfer_2D_V2_OPENMP(double time_step, double time_limit,
                                      double length, double space_step_x,
                                      double width, double space_step_y,
                                      int precision){

    clock_t start_time=clock();
    ll numTimePoint= _cal_num_time(time_step, time_limit);
    ll numSpacePointX= _cal_num_space(length, space_step_x);
    ll numSpacePointY= _cal_num_space(width, space_step_y);
    for (ll t = 0; t < numTimePoint; ++t) {
        for (ll y = 0; y <= numSpacePointY; ++y) {
            for (ll x = 0; x <= numSpacePointX; ++x) {
                _get_value_2D_openmp_v2(time_step, length, space_step_x, width, space_step_y, x, y, t, precision);
            }
        }
    }
    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_OPENMP_V2_without Files is: %f\n",execution_time);
    return 0;

 }