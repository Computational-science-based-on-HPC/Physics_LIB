/**
 * @file thermoutils.h
 * @brief This file contains the utility functions used by the heat transfer simulation in 1D and 2D to help in calculations.
 *
 */

#ifndef PHYSICS_THERMOUTILS_H
#define PHYSICS_THERMOUTILS_H
#define ll long long
#include <stdio.h>

/**
* @brief This is a function calculates the number of steps in the time untils we reach the time limit
* given the rate of change of the time and the time_limit.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @return number of steps we will move as a long long integer.
*/
extern unsigned ll
_cal_num_time(double time_step, double time_limit);

/**
* @brief This is a function calculates the number of steps we will move in the dimentions untils we reach the length
* given the rate of change of the space and the length of the object.
*
* @param length The length of the object.
* @param space_step_x The rate of change of the space in x-axis.
* @return number of steps we will move as a long long integer.
*/
extern ll
_cal_num_space(double length, double space_step);


#endif //PHYSICS_THERMOUTILS_H