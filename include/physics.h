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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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
 * @param gravity gravity of the system
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files will be removed in next update
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

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and each core writes the result to a separate file.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return 0 if there is no error happened interrupted the calculations, each core writes the output to text file named simulate_heat_transfer_1D_MPI_"core num"_ + current time,
* the row represent the time, and the column represent the temperature at this point at that time.
*
* @line
*
* Example:
*  @code{.c}
* #include "physics.h"
*
* int main(void) {
*     heat_equation_1D_P1_MPI(0.01, 0.5, 0.05, 50);
*     return 0;
* }
* @endcode
*/
extern int
heat_equation_1D_P1_MPI(double time_step, double time_limit,
                        double space_step,
                        long long precision);


/**
 * @brief This is a function that simulates the heat transfer in 1D object as wire, and writes the result to a file.
 *
 * In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, we write the result to a file.
 *
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step The rate of change of the space.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, writes the output to text file named simulate_heat_transfer_1D_OPENMP_ + current time,
 * the row represent the time, and the column represent the temperature at this point at that time.
 *
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_1D_P1_OPENMP(0.01, 0.5, 0.05, 50);
 *     return 0;
 * }
 * @endcode
 */
extern int
heat_equation_1D_P1_OPENMP(double time_step, double time_limit,
                           double space_step,
                           long long precision);

////////////////heat equation 2D parallel

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and each core writes the result to a separate file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_MPI_"core num"_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 * 
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_2D_P1_MPI(0.1, 5, 0.1, 0.1, 50);
 *     return 0;
 * }
 * @endcode
 * 
 */

extern int
heat_equation_2D_P1_MPI(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        long long precision);

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and each core writes the result to a separate file.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, we write the result to a file.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return 0 if there is no error happened interrupted the calculations, and write the output to text file named simulate_heat_transfer_2D_OPENMP_ + current time,
 * the output file contains paragraphs each one represent the time slot, each row represent the temperature at this point of the 2D object on x-axis (length), 
 * and the column represent the temperature at this point at that time on y-axis (width).
 * The number of rows in each paragraph (time slot) equals length* space_step_x, and the number of columns equals width* space_step_y.
 * 
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_2D_P1_OPENMP(0.1, 5, 0.1, 0.1, 50);
 *     return 0;
 * }
 * @endcode
 */
extern int
heat_equation_2D_P1_OPENMP(double time_step, double time_limit,
                           double spaceX_step, double spaceY_step,
                           long long precision);

////////////////heat equation 1D serial
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
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_1D_serial(0.01, 0.5, 0.05,50);
 *     return 0;
 * }
 * @endcode
 */

extern int
heat_equation_1D_serial(double time_step, double time_limit, double space_step, long long precision);


////////////////heat equation 2D serial

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
 *
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_2D_serial(0.1, 5, 0.1, 0.1,50);
 *     return 0;
 * }
 * @endcode
 */

extern int
heat_equation_2D_serial(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        long long precision);


////////////////heat equation 1D parallel execution time without i/o

/**
* @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
*
* In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
* Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
* Then, we return the execution time without I/O.
*
* @param time_step The rate of change of the time.
* @param time_limit The time that we want to measure the temperature of the object after.
* @param space_step The rate of change of the space.
* @param precision The number of vectors we use in the calculations.
* @return execution time without I/O.
*
* @line
*
* Example:
*  @code{.c}
* #include "physics.h"
*
* int main(void) {
*     heat_equation_execution_time_1D_P1_MPI(0.01, 0.5, 0.05, 50);
*     return 0;
* }
* @endcode
*/
extern double
heat_equation_execution_time_1D_P1_MPI(double time_step, double time_limit,
                                       double space_step,
                                       long long precision);

/**
 * @brief This is a function that simulates the heat transfer in 1D object as wire, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 1D object as wire and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, return the execution time without I/O.
 *
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step The rate of change of the space.
 * @param precision The number of vectors we use in the calculations.
 * @return execution time without I/O.
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_execution_time_1D_P1_OPENMP(0.01, 0.5, 0.05, 50);
 *     return 0;
 * }
 * @endcode
*/

extern double
heat_equation_execution_time_1D_P1_OPENMP(double time_step, double time_limit,
                                          double space_step,
                                          long long precision);


////////////////heat equation 2D parallel execution time without i/o

/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then we parallelize the function using MPI, as we divide the object into equal parts and each core calculates the temperature of its part.
 * Then, we return the execution time without I/O.
 *
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return execution time without I/O.
 *
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_execution_time_2D_P1_MPI(0.1, 5, 0.1, 0.1,50);
 *     return 0;
 * }
 * @endcode
 */
extern double
heat_equation_execution_time_2D_P1_MPI(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       long long precision);


/**
 * @brief This is a function that simulates the heat transfer in 2D object, and return the execution time without I/O.
 *
 * In this simulation, we simulate heat propagation in 2D object as square or rectangle and the change in its tempreture over time using fourier transform.
 * Then, we parallelize the function using OPENMP, as we provide specific number of threads and make the iterations of the loop are divided into equal-sized chunks, and each chunk is assigned to a thread.
 * Then, return the execution time without I/O.
 * 
 * @param time_step The rate of change of the time.
 * @param time_limit The time that we want to measure the temperature of the object after.
 * @param space_step_x The rate of change of the space in x-axis.
 * @param space_step_y The rate of change of the space in y-axis.
 * @param precision The number of vectors we use in the calculations.
 * @return the execution time without I/O.
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_execution_time_2D_P1_OPENMP(0.1, 5, 0.1, 0.1, 50);
 *     return 0;
 * }
 * @endcode
 * 
 */

extern double heat_equation_execution_time_2D_P1_OPENMP(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        long long precision);


////////////////heat equation 1D serial execution time without i/o

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
 *
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_execution_time_1D_serial(0.01, 0.5, 0.05,50);
 *     return 0;
 * }
 * @endcode
*/

extern double
heat_equation_execution_time_1D_serial(double time_step, double time_limit, double space_step, long long precision);

////////////////heat equation 2D serial execution time without i/o

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
 * 
 * @line
 *
 * Example:
 *  @code{.c}
 * #include "physics.h"
 *
 * int main(void) {
 *     heat_equation_execution_time_2D_serial(0.1, 5, 0.1, 0.1,50);
 *     return 0;
 * }
 * @endcode
 */

extern double
heat_equation_execution_time_2D_serial(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       long long precision);

/**
 * @brief finalize the MPI and de-allocate the resources.
 * @warning This function must be called after using any collection of MPI based functions
 *
 * @line
 *
 * Example:
 * @code{.c}
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