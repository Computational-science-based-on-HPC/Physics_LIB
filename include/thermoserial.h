#ifndef PHYSICS_THERMOSERIAL_H
#define PHYSICS_THERMOSERIAL_H
#include "thermoutils.h"


 extern double
 _get_value_1D(double time_step,
                         double space_step,
                         double x, double t,
                         int precision);

 extern double
 _get_value_2D(double time_step,
               double length, double space_step_x, double width, double space_step_y,
               int x, int y, int t,
               int precision);

 extern int
 _simulate_heat_transfer_1D_serial(double time_step, double time_limit,
                         double length, double space_step,
                         int precision);


 extern int
 _simulate_heat_transfer_2D_serial(double time_step, double time_limit,
                        double length, double space_step_x,
                        double width, double space_step_y,
                        int precision);

extern int
_execution_time_heat_transfer_1D_serial(double time_step, double time_limit,
                                        double length, double space_step,
                                        int precision);
extern int
_execution_time_heat_transfer_2D_serial(double time_step, double time_limit,
                                  double length, double space_step_x,
                                  double width, double space_step_y,
                                  int precision);


#endif //PHYSICS_THERMOSERIAL_H
