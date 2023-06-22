//
// Created by jghal on 6/19/2023.
//

#ifndef PHYSICS_THERMOPARA_H
#define PHYSICS_THERMOPARA_H
#include "thermoutils.h"

extern double
_get_value_1D_mpi(double time_step, 
                 double space_step,
                 double x, double t,
                 int precision);

 extern double
 _get_value_1D_openmp_V1(double time_step,
                         double space_step,
                         double x, double t,
                         int precision);

 extern double
 _get_value_1D_openmp_V2(double time_step,
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

 extern double
 _get_value_2D_openmp_v2(double time_step,
                        double length, double space_step_x,
                        double width, double space_step_y,
                         int x, int y, int t,
                         int precision);

extern int
_simulate_heat_transfer_1D_MPI(double time_step, double time_limit, 
                        double length, double space_step, 
                        int precision);

extern int
_simulate_heat_transfer_1D_OPENMP(double time_step, double time_limit, 
                        double length, double space_step, 
                        int precision);

extern int
_simulate_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit, 
                        double length, double space_step, 
                        int precision);


 extern int
 _simulate_heat_transfer_2D_MPI(double time_step, double time_limit,
                                double length, double space_step_x,
                                double width, double space_step_y,
                                int precision);////////////

 extern int
 _simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
                                   double length, double space_step_x,
                                   double width, double space_step_y,
                                   int precision);

 extern int
 _simulate_heat_transfer_2D_OPENMP_V2(double time_step, double time_limit,
                                      double length, double space_step_x,
                                      double width, double space_step_y,
                                      int precision);



extern int
_execution_time_heat_transfer_1D_MPI(double time_step, double time_limit,
                               double length, double space_step,
                               int precision);

extern int
_execution_time_heat_transfer_1D_OPENMP(double time_step, double time_limit,
                                  double length, double space_step,
                                  int precision);

extern int
_execution_time_heat_transfer_1D_OPENMP_V2(double time_step, double time_limit,
                                     double length, double space_step,
                                     int precision);


extern int
_execution_time_heat_transfer_2D_MPI(double time_step, double time_limit,
                               double length, double space_step_x,
                               double width, double space_step_y,
                               int precision);

//extern int
//_simulate_heat_transfer_2D_OPENMP(double time_step, double time_limit,
//                                  double length, double space_step_x,
//                                  double width, double space_step_y,
//                                  int precision);
//
//extern int
//_simulate_heat_transfer_2D_OPENMP_V2(double time_step, double time_limit,
//                                     double length, double space_step_x,
//                                     double width, double space_step_y,
//                                     int precision);


#endif //PHYSICS_THERMOPARA_H
