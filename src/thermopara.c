/**
 * @file thermopara.c
 * @brief This file contains the implementation of the parallel version of the heat transfer simulation in 1D and 2D.
 *
 * In the 1D heat equation, we make a model that thermal energy flow in a 1D object so we place the object along the x-axis.
 * We show that an object has a fixed length and lies along an interval [0, 1] under the condition of both time and thermal energy we have some conditions that affect the temperature that is time and point of x on the object.
 * We also have thermal diffusivity that shows the rate of the temperature to be transferred to the object(how quickly the heat moves throw the object)
 * As, we separated variables of the thermal energy problem to calculate u(x, t) we make this function separated to the function of x-axis X(x) and function of time T(t)
 * u(x, t) = X(t) T(t)
 * Then, we calculate the change throw the time t to get the equation
 * u(x,t)= ∑b(n) * sin(nΠx/L)
 * b(n)= (2/L) ∫f(x)*sin(nΠx/L)dx
 * In addition, parallelizing this functions will be done by using OpenMP and MPI.
 * 
 * In the 2D heat equation,  we make a model that thermal energy flow in a 2D object, so we place the object on xy-plane under the condition of time.
 * We represent the object by [0, a] x [0, b] as a is the length on x-axis and b is the width on y-axis, so, the tempertaure affected by function of x-axis X(x) and and y-axis Y(y) and time T(t).
 * represented by u(x, y, t) = X(x) Y(y) T(t)
 * Then, we separated the variables as 1D to get in each one of them the constant change as we get u(t) by using the diffusivity on u(xx) and u(yy)
 * u(t) = c^2 (u(xx) + u(yy))
 * So, we get final equation to be as
 * u(x,y,t)= ∑∑A(mn) sin(μ(m)x)sin(υ(n)y) exp(-λ(mn)t)
 * μ(m)= mΠ/a , υ(n)=nΠ/b , λ(mn) = c√μ(m)^2 + υ(n)^2 
 * A(mn)= 4/ab∫∫f(x,y)sin(μ(m)x)sin(υ(n)y) dy dx
 * In addition, parallelizing this functions will be done by using OpenMP and MPI.
 */

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

/**
* @brief This is a function calculates the value of specific point in the space at specific time in 1D.
*
* In this function, we calculate the temperature of point x at time t given the rate of change of x and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
*
* @param time_step The rate of change of the time.
* @param space_step The rate of change of the space.
* @param x The point we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/
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

/**
* @brief This is a function calculates the value of specific point in the space at specific time in 1D.
*
* In this function, we calculate the temperature of point x at time t given the rate of change of x and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
* Then we parallelize the for loop using OpenMP as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
* the shared data that all the threads need to access is the sum, x_real, t_real, precision.
* the private data that each thread need to access is the exponential, spaceXTerm, coeff.
*
* @param time_step The rate of change of the time.
* @param space_step The rate of change of the space.
* @param x The point we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/
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

/**
* @brief This is a function calculates the value of specific point in the space at specific time in 1D.
*
* In this function, we calculate the temperature of point x at time t given the rate of change of x and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
* Then we parallelize the for loop using OpenMP as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
* Each section block represents a distinct section of code that will be executed by a separate thread. In this case, we have three sections
*
* @param time_step The rate of change of the time.
* @param space_step The rate of change of the space.
* @param x The point we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/
double _get_value_1D_openmp_V2(double time_step, double space_step, double x, double t, int precision){
    double sum = 0.0, exponential, spaceXTerm, coeff;
    double x_real = x * space_step;
    double t_real = t * time_step;

    #pragma omp parallel num_threads(THREADS)
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

/**
* @brief This is a function calculates the value of specific point in the space at specific time in 2D.
*
* In this function, we calculate the temperature of point x and y at time t given the rate of change of x, y and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
*
* @param time_step The rate of change of the time.
* @param length The length of the object.
* @param space_step_x The rate of change of the space in x-axis.
* @param width The width of the object.
* @param space_step_y The rate of change of the space in y-axis.
* @param x The point in x-axis we want to calculate the temperature at.
* @param y The point in y-axis we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/
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

/**
* @brief This is a function calculates the value of specific point in the space at specific time in 2D.
*
* In this function, we calculate the temperature of point x and y at time t given the rate of change of x, y and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
* Then we parallelize the for loop using OpenMP as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
* the shared data that all the threads need to access is the sum, x_real, t_real, precision.
* the private data that each thread need to access is the exponential, spaceXTerm, coeff.
*
* @param time_step The rate of change of the time.
* @param length The length of the object.
* @param space_step_x The rate of change of the space in x-axis.
* @param width The width of the object.
* @param space_step_y The rate of change of the space in y-axis.
* @param x The point in x-axis we want to calculate the temperature at.
* @param y The point in y-axis we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/
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


/**
* @brief This is a function calculates the value of specific point in the space at specific time in 2D.
*
* In this function, we calculate the temperature of point x and y at time t given the rate of change of x, y and t, using fourier transform
* and the number of times we sum the series depend on the number of precision provided.
* Then we parallelize the for loop using OpenMP as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
* Each section block represents a distinct section of code that will be executed by a separate thread. In this case, we have three sections
*
* @param time_step The rate of change of the time.
* @param length The length of the object.
* @param space_step_x The rate of change of the space in x-axis.
* @param width The width of the object.
* @param space_step_y The rate of change of the space in y-axis.
* @param x The point in x-axis we want to calculate the temperature at.
* @param y The point in y-axis we want to calculate the temperature at.
* @param t The current time we are at.
* @param precision The number of vectors we use in the calculations.
* @return sum we calculated.
*/

double _get_value_2D_openmp_v2(double time_step,
                        double length, double space_step_x, double width, double space_step_y,
                        int x, int y, int t,
                        int precision){
    double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
    double x_real = x * space_step_x;
    double y_real = y * space_step_y;
    double t_real = t * time_step;

    #pragma omp parallel num_threads(THREADS)
    for (ll m = 1; m < precision; ++m)
    {
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

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and each core writes the result to a separate file.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return 0 if there is no error happened interrupted the calculations, each core writes the output to text file named simulate_heat_transfer_1D_MPI_"core num"_ + current time,
* the row represent the time, and the column represent the temperature at this point at that time.
*/

int _simulate_heat_transfer_1D_MPI(double time_step, double time_limit, double space_step, int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
    double length =10.0;
     MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_1D_MPI_1_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr1 = fopen(_file_name, "w");
    sprintf(_file_name, "simulate_heat_transfer_1D_MPI_2_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr2 = fopen(_file_name, "w");
    sprintf(_file_name, "simulate_heat_transfer_1D_MPI_3_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr3 = fopen(_file_name, "w");
//    fptr1 = fopen("1_1D_MPI.txt", "w");
//    fptr2 = fopen("2_1D_MPI.txt", "w");
//    fptr3 = fopen("3_1D_MPI.txt", "w");
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

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and writes the result to a file.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
* Then, we write the result to a file.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return 0 if there is no error happened interrupted the calculations, writes the output to text file named simulate_heat_transfer_1D_OPENMP_ + current time,
* the row represent the time, and the column represent the temperature at this point at that time.
*/

int _simulate_heat_transfer_1D_OPENMP(double time_step, double time_limit, double space_step, int precision){
    clock_t start_time=clock();
    double length =10.0;
    FILE *fptr;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_1D_OPENMP_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr = fopen(_file_name, "w");

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

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and writes the result to a file.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using OPENMP,as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
* and then, write the result to a file.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return 0 if there is no error happened interrupted the calculations, writes the output to text file named simulate_heat_transfer_1D_OPENMP_V2_ + current time,
* the row represent the time, and the column represent the temperature at this point at that time.
*/
int
_simulate_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit, double space_step, int precision){
    clock_t start_time=clock();
    double length =10.0;
    FILE *fptr;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_1D_OPENMP_V2_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr = fopen(_file_name, "w");


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

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and each core writes the result to a separate file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_MPI_"core num"_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 */
 int
 _simulate_heat_transfer_2D_MPI(double time_step, double time_limit,
                                double space_step_x,
                                double space_step_y,
                                int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
     double length =2.0, width =2.0;
     MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(NULL, NULL);
    }
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    char _file_name[2076];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    sprintf(_file_name, "simulate_heat_transfer_2D_MPI_1_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr1 = fopen(_file_name, "w");
    sprintf(_file_name, "simulate_heat_transfer_1D_MPI_2_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr2 = fopen(_file_name, "w");
    sprintf(_file_name, "simulate_heat_transfer_1D_MPI_3_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fptr3 = fopen(_file_name, "w");
//    fptr1 = fopen("1_2D_MPI.txt", "w");
//    fptr2 = fopen("2_2D_MPI.txt", "w");
//    fptr3 = fopen("3_2D_MPI.txt", "w");

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
            }
            if(my_rank == 1){
                fprintf(fptr1, "\n\n");
            }else if(my_rank == 2){
                fprintf(fptr2, "\n\n");
            }else if(my_rank == 3){
                fprintf(fptr3, "\n\n");
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
        printf("The time taken in MPI_2d With I/O is: %f\n", end_time - start_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
            
    return 0;
 }


/**
 * @brief This is a function that simulates the heat transfer in 2D object, and each core writes the result to a separate file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, we write the result to a file.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_OPENMP_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 */

 int
 _simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double space_step_x,
                                   double space_step_y,
                                   int precision){

    clock_t start_time=clock();
    double length =2.0, width =2.0;
    FILE *fptr;
     char _file_name[2076];
     time_t t = time(NULL);
     struct tm tm = *localtime(&t);
     sprintf(_file_name, "simulate_heat_transfer_2D_OPENMP_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
     fptr = fopen(_file_name, "w");

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

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and each core writes the result to a separate file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using OPENMP,as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
 * and then, write the result to a file.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_OPENMP_V2_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 */

 int
 _simulate_heat_transfer_2D_OPENMP_V2(double time_step, double time_limit,
                                      double space_step_x,
                                      double space_step_y,
                                      int precision){

    clock_t start_time=clock();
    double length =2.0, width =2.0;
    FILE *fptr;
     char _file_name[2076];
     time_t t = time(NULL);
     struct tm tm = *localtime(&t);
     sprintf(_file_name, "simulate_heat_transfer_2D_OPENMP_V2_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
     fptr = fopen(_file_name, "w");



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

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
* Then, we return the execution time without I/O.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return execution time without I/O.
*/

double _execution_time_heat_transfer_1D_MPI(double time_step, double time_limit,
                                     double space_step,
                                     int precision){
    // MPI_Init(NULL, NULL);
    int initialized, finalized;
    double length =10.0;
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
    double execution_time;
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;
        printf("The time taken in MPI_1d Without I/O is: %f\n", execution_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
        
    return execution_time;
}

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
* Then, return the execution time without I/O.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return execution time without I/O.
*/

double _execution_time_heat_transfer_1D_OPENMP(double time_step, double time_limit,
                                        double space_step,
                                        int precision){

    clock_t start_time=clock();
    double length =10.0;
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
    return execution_time;
}

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using OPENMP,as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
* and then, return the execution time without I/O.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return execution time without I/O.
*/
double _execution_time_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit,
                                           double space_step,
                                           int precision){

    clock_t start_time=clock();
    double length =10.0;
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
    return execution_time;

}

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return execution time without I/O.
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 */

double _execution_time_heat_transfer_2D_MPI(double time_step, double time_limit,
                                     double space_step_x,
                                     double space_step_y,
                                     int precision){
    
    int initialized, finalized;
    double length =2.0, width =2.0;
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
    double execution_time;
    if (my_rank == 0)
    {
        end_time = MPI_Wtime();
        execution_time = end_time - start_time;
        printf("The time taken in MPI_2d Without I/O is: %f\n", execution_time);
    }
    // MPI_Finalize();
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
    
    return execution_time;

}

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, return the execution time without I/O.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return the execution time without I/O.
 */

 double
 _execution_time_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double space_step_x,
                                   double space_step_y,
                                   int precision){

    clock_t start_time=clock();
     double length =2.0, width =2.0;
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
    return execution_time;

 }

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using OPENMP,as we provide specific number of threads, then we used the parallel sections which creates the parallel regions where the sections can be executed in parallel. 
 * and then, write the result to a file.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return the execution time without I/O.
 */

 double
 _execution_time_heat_transfer_2D_V2_OPENMP(double time_step, double time_limit,
                                      double space_step_x,
                                      double space_step_y,
                                      int precision){

    clock_t start_time=clock();
     double length =2.0, width =2.0;
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
    return execution_time;

 }