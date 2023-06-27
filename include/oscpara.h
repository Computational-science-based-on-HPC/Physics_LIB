//
// Created by jghal on 6/16/2023.
//
/**
 * @file oscpara.h
 * @brief This file contains the implementation of the parallel versions of the oscillation simulation in 1D and 2D.
 *
 */
#ifndef PHYSICS_OSCPARA_H
#define PHYSICS_OSCPARA_H
/**
 *  @brief function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
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
extern int
_simulate_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);
/**
 *  @brief This function calculate the execution time of simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
 *  @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
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
 * @return
 */
extern int
_execution_time_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                           double Ao,
                                           double Vo, double FI,
                                           double time_limit, double step_size, double damping_coefficent,
                                           int number_of_files);
/**
 *  @brief This function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
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
 * @return
 */
extern int
_simulate_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                 double Vo, double FI,
                                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  @brief function calculate execution time simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
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
 * @return
 */
extern int
_execution_time_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k,
                                       double Ao,
                                       double Vo, double FI,
                                       double time_limit, double step_size, double damping_coefficent,
                                       int number_of_files);
/**
 *  @brief This function calculate execution time simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using Openmp, The for loop is divided upon multiple threads.
 *  @note tried using sections pragma ended up applying more overhead in the code so performance decreased
 * @warning don't set gravity less than 0.
 * @warning if max_amplitude is greater than the length the max_amplitude automatically set to the length value
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
 * @param num_of_threads number of threads needed to execute the code
 * @return
 */
extern double
_execution_time_damped_os_parallel_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                          double Ao,
                                          double Vo, double FI,
                                          double time_limit, double step_size, double damping_coefficent,
                                          int number_of_files, int num_of_threads);

#endif // PHYSICS_OSCPARA_H
