#include "../include/thermoserial.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/thermoutils.h"

#define  ll long long
// #define M_PI 3.14159265358979323846264338327


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