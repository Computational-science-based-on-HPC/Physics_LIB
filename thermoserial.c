#include "thermoserial.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "thermoutils.h"

#define  ll long long
#define M_PI 3.14159265358979323846264338327

// double _get_value_1D(struct TimeParam* time, struct SpaceParam* space, double x, double t, int precision)
// {
//     double sum = 0.0, exponential, spaceXTerm, coeff;
//     double x_real = x * space->delta_x;
//     double t_real = t * time->delta_t;

//     for (int k = 0; k < precision; k++)
//     {
//         exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
//         spaceXTerm = sin((double)(2 * k + 1) * M_PI * x_real / 2);
//         coeff = 1 / (2 * k + 1);
//         sum += coeff * exponential * spaceXTerm;
//     }

//     sum *= 200 / M_PI;
//     return sum;
// }

// double _get_value_2D(struct TimeParam* time, struct SpaceParam2D* space, int x, int y, int t, int precision) {
//     double sum = 0.0, exponential, spaceXTerm, spaceYTerm, coeff;
//     double x_real = x * space->delta_x;
//     double y_real = y * space->delta_y;
//     double t_real = t * time->delta_t;
//     for (ll m = 1; m < precision; ++m) {
//         for (ll n = 1; n < precision; ++n) {
//             exponential = exp(-(M_PI * M_PI) * (m * m + n * n) * t_real / 36);
//             spaceXTerm = sin((double) m * M_PI * x_real / space->length);
//             spaceYTerm = sin((double) n * M_PI * y_real / space->width);
//             // Find Amn constant and multiply it with the sum
//             coeff = (1 + pow(-1, m + 1)) * (1 - cos(n * M_PI / 2)) / (m * n);
//             sum += coeff * exponential * spaceXTerm * spaceYTerm;
//         }
//     }
//     sum *= 200 / (M_PI * M_PI);
//     return sum;
// }

// int _simulate_heat_transfer_1D_serial(struct TimeParam* time_param, struct SpaceParam* space_param,  int precision){
//     FILE *fptr;
//     fptr = fopen("1D_serial_V1.txt", "w");

//     ll numTimePoint;
//     _cal_num_time(time_param, &numTimePoint);

//     ll numSpacePoint;
//     _cal_num_space(space_param, &numSpacePoint);

//     for (int t = 0; t < numTimePoint; t++) {
//         for (int x = 0; x < numSpacePoint; x++) {
//             fprintf(fptr, "%f ", _get_value_1D(time_param, space_param, x, t, precision));
//         }
//         fprintf(fptr, "\n");
//     }

//     fclose(fptr);
//     return 0;

// }

// int _simulate_heat_transfer_2D_serial(struct TimeParam* time_param, struct SpaceParam2D* space_param, struct TempParam* temp_param, int precision){
//     FILE *fptr;
//     fptr = fopen("2D_serial_V1.txt", "w");
    
//     ll numTimePoint;
//     ll numSpacePointX;
//     ll numSpacePointY;

//     _cal_num_time(time_param, &numTimePoint);
//     _cal_num_space_2D(time_param, space_param, &numSpacePointX, &numSpacePointY);
    
//     for (ll t = 0; t < numTimePoint; ++t) {
//         for (ll y = 1; y < numSpacePointY; ++y) {
//             for (ll x = 1; x < numSpacePointX; ++x) {
//                 fprintf(fptr, "%f ", _get_value_2D(time_param, space_param, x, y, t, precision));
//             }
//             fprintf(fptr, "\n");
//         }
//         fprintf(fptr, "\n\n");
//     }

//     fclose(fptr);
//     return 0;
// }
