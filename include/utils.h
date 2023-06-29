//
// Created by jghal on 6/16/2023.
//
/**
 * @file utils.h
 * @brief This file contains utility functions used by oscillations simulation to help in calculations.
 *
 */
#ifndef PHYSICS_UTILS_H
#define PHYSICS_UTILS_H

/**
 * @brief Validate on input entered to oscillation.
 *
 * @param x distance x
 * @param y distance y
 * @param length max length
 * @param mass  mass of bob
 * @param gravity
 * @param k stiffeness of spring
 * @param time_limit how much time needed to simulate
 * @param step_size how much time change per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @param Fo
 * @return 1 if all are valid, 0 if [mass < 0 || k < 0 || time_limit <= 0 || step_size < 0 || length <= 0 || gravity <= 0 || damping_coefficent < 0], and -1 if max_length > length
 */
int _valid_osc(double x, double y, double length, double mass, double gravity, double k, double time_limit,
               double step_size,
               double damping_coefficent, int number_of_files, double Fo);
/**
 * @brief get the minimum between x,y
 *
 * @param x
 * @param y
 * @return min(x,y)
 */
int _min_int(int x, int y);
/**
 * @brief round x value to the nearest int
 * @param x
 * @return rounded x
 */
int _round(double x);
/**
 * @brief calculate dx for elastic pendulum system
 * @param dx
 * @return
 */
double
_dx(double dx);
/**
 * @brief calculate dy for elastic pendulum system
 * @param dy
 * @return
 */
double
_dy(double dy);
/**
 * @brief calculate first function to solve ode in elastic pendulum
 * @param x
 * @param y
 * @param dx
 * @param tx
 * @param ty
 * @param k
 * @param m
 * @param b
 * @param r
 * @return
 */
double
_f1(double x, double y, double dx, double tx, double ty, double k, double m, double b, double r);
/**
 * @brief calculate second function to solve ode in elastic pendulum
 * @param x
 * @param y
 * @param dy
 * @param tx
 * @param ty
 * @param k
 * @param m
 * @param b
 * @param r
 * @param g
 * @return
 */
double
_f2(double x, double y, double dy, double tx, double ty, double k, double m, double b, double r, double g);
void printmemsizestream(char *str, unsigned long ramsize);
/**
 * @brief this function prints memory details for the computer
 * 
 * @return int 
 */
int printmemstream();
/**
 * @brief this function prints the cpu info for the computer
 * 
 * @return int 
 */
int cpu_inf_stream();

#endif // PHYSICS_UTILS_H
