#include "../include/thermoserial.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/thermoutils.h"

#define  ll long long
// #define M_PI 3.14159265358979323846264338327


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

double _get_value_1D(double time_step,
                    double space_step,
                    double x, double t,
                    int precision)
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
double _get_value_2D(double time_step,
                    double length, double space_step_x, double width, double space_step_y,
                    int x, int y, int t,
                    int precision) {
    double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
    double x_real = x * space_step_x;
    double y_real = y * space_step_y;
    double t_real = t * time_step;
    for (ll m = 1; m < precision; ++m) {
        for (ll n = 1; n < precision; ++n) {
            exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
            spaceXTerm = sin((double) m * M_PI * x_real / length);
            spaceYTerm = sin((double) n * M_PI * y_real / width);
            // Find Amn constant and multiply it with the sum
            coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
            sum += coeff * exponential * spaceXTerm * spaceYTerm;
        }
    }
    sum *= 200 / (M_PI * M_PI);
    return sum;
}
/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and write the result to a file.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_1D_serial_ + current time,
* the row represent the time, and the column represent the temperature at this point at that time.
*/

 int _simulate_heat_transfer_1D_serial(double time_step, double time_limit,
                                       double space_step,
                                       int precision){
     clock_t start_time=clock();
     double length =10.0;
     FILE *fptr;
     char _file_name[2076];
     time_t t = time(NULL);
     struct tm tm = *localtime(&t);
     sprintf(_file_name, "simulate_heat_transfer_1D_serial_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
     fptr = fopen(_file_name, "w");

     ll numTimePoint= _cal_num_time(time_step, time_limit);

     ll numSpacePoint= _cal_num_space(length, space_step);

     for (int t = 0; t < numTimePoint; t++) {
         for (int x = 0; x <= numSpacePoint; x++) {
             fprintf(fptr, "%f ", _get_value_1D(time_step, space_step, x, t, precision));
         }
         fprintf(fptr, "\n");
     }

     fclose(fptr);
     clock_t end_time=clock();
     double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
     printf("The value of execution_time 1D_serial_simulation_withFiles is: %f\n",execution_time);
     return 0;

 }

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and write the result to a file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 *
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_serial_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 */

 int _simulate_heat_transfer_2D_serial(double time_step, double time_limit,
                                       double space_step_x,
                                       double space_step_y,
                                       int precision){
     clock_t start_time=clock();
     double length =2.0, width =2.0;
     FILE *fptr;
     char _file_name[2076];
     time_t t = time(NULL);
     struct tm tm = *localtime(&t);
     sprintf(_file_name, "simulate_heat_transfer_2D_serial_%d-%02d-%02d %02d:%02d:%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
     fptr = fopen(_file_name, "w");


     ll numTimePoint= _cal_num_time(time_step, time_limit);
     ll numSpacePointX= _cal_num_space(length, space_step_x);
     ll numSpacePointY= _cal_num_space(width, space_step_y);

     for (ll t = 0; t < numTimePoint; ++t) {
         for (ll y = 0; y <= numSpacePointY; ++y) {
             for (ll x = 0; x <= numSpacePointX; ++x) {
                 fprintf(fptr, "%f ", _get_value_2D(time_step,  length, space_step_x, width, space_step_y, x, y, t, precision));
             }
             fprintf(fptr, "\n");
         }
         fprintf(fptr, "\n\n");
     }

     fclose(fptr);
     clock_t end_time=clock();
     double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
     printf("The value of execution_time 2D_serial_simulation_withFiles is: %f\n",execution_time);
     return 0;
 }


/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return execution time without I/O if there is no error happened interrupted the calculations.
*/

double _execution_time_heat_transfer_1D_serial(double time_step, double time_limit,
                                            double space_step,
                                            int precision){
    clock_t start_time=clock();
    double length =10.0;
    ll numTimePoint= _cal_num_time(time_step, time_limit);

    ll numSpacePoint= _cal_num_space(length, space_step);

    for (int t = 0; t < numTimePoint; t++) {
        for (int x = 0; x <= numSpacePoint; x++) {
            _get_value_1D(time_step, space_step, x, t, precision);
        }
    }

    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 1D_serial is: %f\n",execution_time);
    return execution_time;

}


/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 *
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return execution time without I/O if there is no error happened interrupted the calculations.
 */
double _execution_time_heat_transfer_2D_serial(double time_step, double time_limit,
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
                 _get_value_2D(time_step, length, space_step_x, width, space_step_y, x, y, t, precision);
            }
        }
    }

    clock_t end_time=clock();
    double execution_time=(double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("The value of execution_time 2D_serial is: %f\n",execution_time);
    return execution_time;


}