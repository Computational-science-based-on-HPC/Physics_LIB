/**
 * @file thermoserial.h
 * @brief This file contains the implementation of the serial version of the heat transfer simulation in 1D and 2D.
 *
 * In the 1D heat equation, we make a model that thermal energy flow in a 1D object so we place the object along the x-axis.
 * We show that an object has a fixed length and lies along an interval [0, 1] under the condition of both time and thermal energy we have some conditions that affect the temperature that is time and point of x on the object.
 * We also have thermal diffusivity that shows the rate of the temperature to be transferred to the object(how quickly the heat moves throw the object)
 * As, we separated variables of the thermal energy problem to calculate u(x, t) we make this function separated to the function of x-axis X(x) and function of time T(t)
 * u(x, t) = X(t) T(t)
 * Then, we calculate the change throw the time t to get the equation
 * u(x,t)= ∑b(n) * sin(nΠx/L)
 * b(n)= (2/L) ∫f(x)*sin(nΠx/L)dx
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
 */


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
                         double space_step,
                         int precision);


 extern int
 _simulate_heat_transfer_2D_serial(double time_step, double time_limit,
                        double space_step_x,
                        double space_step_y,
                        int precision);

extern double
_execution_time_heat_transfer_1D_serial(double time_step, double time_limit,
                                        double space_step,
                                        int precision);
extern double
_execution_time_heat_transfer_2D_serial(double time_step, double time_limit,
                                  double space_step_x,
                                  double space_step_y,
                                  int precision);


#endif //PHYSICS_THERMOSERIAL_H
