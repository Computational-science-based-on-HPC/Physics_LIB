#ifndef PHYSICS_THERMOPARA_H
#define PHYSICS_THERMOPARA_H
#include "thermoutils.h"

extern double
_get_value_1D_mpi(double time_step, 
                 double space_step,
                 double x, double t,
                 int precision);

 extern double
 _get_value_1D_openmp(double time_step,
                         double space_step,
                         double x, double t,
                         int precision);

 extern double
 _get_value_2D_mpi(double time_step,
                        double length, double space_step_x,
                        double width, double space_step_y,
                         int x, int y, int t,
                         int precision);

 extern double
 _get_value_2D_openmp(double time_step,
                        double length, double space_step_x,
                        double width, double space_step_y,
                         int x, int y, int t,
                         int precision);

extern int
_simulate_heat_transfer_1D_MPI(double time_step, double time_limit, 
                        double space_step,
                        int precision);

extern int
_simulate_heat_transfer_1D_OPENMP(double time_step, double time_limit, 
                        double space_step,
                        int precision);

 extern int
 _simulate_heat_transfer_2D_MPI(double time_step, double time_limit,
                                double space_step_x,
                                double space_step_y,
                                int precision);

 extern int
 _simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double space_step_x,
                                   double space_step_y,
                                   int precision);

extern double
_execution_time_heat_transfer_1D_MPI(double time_step, double time_limit,
                               double space_step,
                               int precision);

extern double
_execution_time_heat_transfer_1D_OPENMP(double time_step, double time_limit,
                                  double space_step,
                                  int precision);


extern double
_execution_time_heat_transfer_2D_MPI(double time_step, double time_limit,
                               double space_step_x,
                               double space_step_y,
                               int precision);

extern double
_execution_time_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double space_step_x,
                                   double space_step_y,
                                   int precision);


#endif //PHYSICS_THERMOPARA_H
