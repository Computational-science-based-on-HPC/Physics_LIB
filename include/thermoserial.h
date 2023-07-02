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

 extern long double
 _get_value_1D(double time_step,
                         double space_step,
                         long long x, unsigned long long t,
                         long long precision);

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
 extern long double
 _get_value_2D(double time_step,
               double length, double space_step_x, double width, double space_step_y,
               long long x, long long y, unsigned long long t,
               long long precision);

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

 extern int
 _simulate_heat_transfer_1D_serial(double time_step, double time_limit,
                         double space_step,
                         long long precision);



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
 extern int
 _simulate_heat_transfer_2D_serial(double time_step, double time_limit,
                        double space_step_x,
                        double space_step_y,
                        long long precision);

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

extern double
_execution_time_heat_transfer_1D_serial(double time_step, double time_limit,
                                        double space_step,
                                        long long precision);

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

extern double
_execution_time_heat_transfer_2D_serial(double time_step, double time_limit,
                                  double space_step_x,
                                  double space_step_y,
                                  long long precision);


#endif //PHYSICS_THERMOSERIAL_H
