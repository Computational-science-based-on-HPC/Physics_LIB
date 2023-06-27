/**
 * @file physics.h
 * @brief This file contains collection of all simulations calls.
 *
 */
#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H
/**
 *  @brief function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm
 *
 * Output of this simulation will be resulted in file named "damped_os_serial_displacement_YYYY-MM-DD HH:MM:SS.txt"
 *
 *
 *
 *
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
 * @param number_of_files currently nulled
 * @return
 *
 * @line
 *
 * Example:
 * @code{.c}
 * #include "physics.h"
 *
 *  int main(void) {
 *     damped_os_serial(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */
extern int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                 double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  @brief function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
 *
 * Output of this simulation will be resulted in file named "PROCESS-RANK_damped_os_parallel_v1_displacement_YYYY-MM-DD HH:MM:SS.txt"
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
 * @return
 * @line
 *
 * Example:
 * @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_v1(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */
extern int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  @brief This function simulates simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
 *
 * Output of this simulation will be resulted in file named "PROCESS-RANK_damped_os_parallel_v2_displacement_YYYY-MM-DD HH:MM:SS.txt"
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
 * @param number_of_files currently nulled
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_v2(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */
extern int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 * @brief function simulates the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 *
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
 *
 * Output of this simulation will be resulted in file named "elastic_pendulum_x_YYYY-MM-DD HH:MM:SS.txt" and "elastic_pendulum_y_YYYY-MM-DD HH:MM:SS.txt"
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
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     elastic_pendulum(2.0, 5.0, 0.5, 9.81, 0.005, 0, 2.0, 0.0, 0, 100000, 0.01, 0, 3)
 *     return 0;
 * }
 * @endcode
 */
extern int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  @brief This function calculate the execution time of simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI and Openmp together, The Number of iteration are divided upon number of processes using MPI, while each processes is calculating its values Openmp used to run the calculations per process in parallel way
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
 * @param number_of_files currently nulled
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_execution_time_v1(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */
extern double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);
/**
 *  @brief This function calculate execution time simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using MPI, The Number of iteration are divided upon number of processes using MPI.
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
 * @param number_of_files currently nulled
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_execution_time_v2(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */

extern double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
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
 *
 *  @note tried using sections pragma ended up applying more overhead in the code so performance decreased
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
 * @param number_of_files currently nulled
 * @param num_of_threads number of threads needed to execute the code
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_execution_time_v3(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 *
 */
extern double
damped_os_parallel_execution_time_v3(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files,int num_of_threads);
/**
 *  @brief This function calculates execution time of simulating simple harmonic motion (Simple Spring Motion).
 *
 *  using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm
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
 * @param number_of_files currently nulled
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_serial_execution(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     return 0;
 * }
 * @endcode
 */
extern double
damped_os_serial_execution(double max_amplitude, double length, double mass, double gravity, double k,
                           double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent,
                           int number_of_files);
/**
 * @brief This function calculates the execution time of simulating the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 *
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
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
 * @return
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     elastic_pendulum_execution(2.0, 5.0, 0.5, 9.81, 0.005, 0, 2.0, 0.0, 0, 100000, 0.01, 0, 3)
 *     return 0;
 * }
 * @endcode
 */
extern double
elastic_pendulum_execution(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);


////////////////heat equation 1D parallel
extern int
heat_equation_1D_P1_MPI(double time_step, double time_limit,
                        double space_step,
                        int precision);

extern int
heat_equation_1D_P1_OPENMP(double time_step, double time_limit,
                           double space_step,
                           int precision);

////////////////heat equation 2D parallel
extern int
heat_equation_2D_P1_MPI(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);

extern int
heat_equation_2D_P1_OPENMP(double time_step, double time_limit,
                           double spaceX_step, double spaceY_step,
                           int precision);

////////////////heat equation 1D serial
extern int
heat_equation_1D_serial(double time_step, double time_limit, double space_step, int precision);


////////////////heat equation 2D serial
extern int
heat_equation_2D_serial(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);


////////////////heat equation 1D parallel execution time without i/o
extern double
heat_equation_execution_time_1D_P1_MPI(double time_step, double time_limit,
                                       double space_step,
                                       int precision);

extern double
heat_equation_execution_time_1D_P1_OPENMP(double time_step, double time_limit,
                                          double space_step,
                                          int precision);


////////////////heat equation 2D parallel execution time without i/o
extern double
heat_equation_execution_time_2D_P1_MPI(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       int precision);

extern double heat_equation_execution_time_2D_P1_OPENMP(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);


////////////////heat equation 1D serial execution time without i/o
extern double
heat_equation_execution_time_1D_serial(double time_step, double time_limit, double space_step, int precision);

////////////////heat equation 2D serial execution time without i/o
extern double
heat_equation_execution_time_2D_serial(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       int precision);

/**
 * @brief finalize the MPI and de-allocate the resources.
 * @warning This function must be called after using any collection of MPI based functions
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     damped_os_parallel_execution_time_v2(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     damped_os_parallel_execution_time_v1(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     damped_os_parallel_execution_time_v3(10.0, 14.0, 1.0, 9.8, 1.0, -1.0, 0.0,-0.1, 40000000, 0.01, 0.1, 3);
 *     finalize();
 *     return 0;
 * }
 * @endcode
 */
extern void
finalize();

#endif // PHYSICS_PHYSICS_H