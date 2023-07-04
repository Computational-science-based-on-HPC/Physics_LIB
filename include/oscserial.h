//
// Created by jghal on 6/16/2023.
//

/**
 * @file oscserial.h
 * @brief This file contains the implementation of the serial version of the oscillation simulation in 1D and 2D.
 *
 */
#ifndef PHYSICS_OSCSERIAL_H
#define PHYSICS_OSCSERIAL_H
/**
 *  @brief This function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm.
 *
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return integer value 0 if the program executes without any errors -1 if there is an error occurred during calculations
 */
extern char*
_simulate_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  @brief This function calculates execution time of simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm.
 *
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return execution time of simulation
 */
extern double
_execution_time_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 * @brief This function simulates the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 *
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
 *
 * @warning don't set gravity less than 0.
 * @warning if sqrt(Yo^2 + Xo^2) is greater than the length the Xo = sqrt(length ^ 2 - Yo ^ 2).
 *
 * @param r rest length of spring
 * @param length max length of spring
 * @param mass mass of bob suspended in spring
 * @param gravity
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return integer value 0 if the program executes without any errors -1 if there is an error occurred during calculations
 */
extern char*
_simulate_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 * @brief This function calculates the execution time of simulating the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 *
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
 *
 * @warning don't set gravity less than 0.
 * @warning if sqrt(Yo^2 + Xo^2) is greater than the length the Xo = sqrt(length ^ 2 - Yo ^ 2).
 *
 * @param r rest length of spring
 * @param length max length of spring
 * @param mass mass of bob suspended in spring
 * @param gravity
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return execution time if simulation
 */
extern double
_execution_time_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);

#endif //PHYSICS_OSCSERIAL_H
