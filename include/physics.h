#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H
/**
 *  This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                 double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 * This function simulates the motion of (elastic pendulum/2D-spring/spring pendulum) system.
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
 */
extern int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
/**
 *  This function calculate the execution time of simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);
/**
 *  This function calculate execution time simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);
/**
 *  This function calculate execution time simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using Openmp, The for loop is divided upon multiple threads.
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
 */
extern double
damped_os_parallel_execution_time_v3(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files,int num_of_threads);
/**
 *  This function calculates execution time of simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
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
 */
extern double
damped_os_serial_execution(double max_amplitude, double length, double mass, double gravity, double k,
                           double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent,
                           int number_of_files);
/**
 * This function calculates the execution time of simulating the motion of (elastic pendulum/2D-spring/spring pendulum) system.
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
extern void 
finalize();

#endif // PHYSICS_PHYSICS_H