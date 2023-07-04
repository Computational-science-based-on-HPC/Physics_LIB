/**
 * @file physics.c
 * @brief This file contains collection of all simulations calls.
 *
 */
#include "../include/physics.h"
#include "../include/oscserial.h"
#include "../include/oscpara.h"
#include "../include/thermopara.h"
#include "../include/thermoserial.h"
#include "../include/utils.h"
#include <mpi.h>
/**
 *  This function used to call _simulate_damped_os_serial
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
char *damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                       double FI,
                       double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Damped oscillation serial");
    return _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                      time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 *  This function used to call _simulate_damped_os_parallel_mpi_omp
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
char *damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                            double FI,
                            double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Damped oscillation mpi 1");
    return _simulate_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                                time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 * This function used to call _simulate_damped_os_parallel_mpi
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
char *damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                            double FI,
                            double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Damped oscillation mpi 2");
    return _simulate_damped_os_parallel_mpi(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                            time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 * This function used to call _simulate_elastic_pendulum
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
char *elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                       double Yo,
                       double Vo,
                       double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Elastic pendulum");
    return _simulate_elastic_pendulum(r, length, mass, gravity, k, Ao, Xo,
                                      Yo,
                                      Vo,
                                      time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 * This function used to call _execution_time_damped_os_parallel_mpi_omp
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
double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
    _mkdir("Logs");
    _mkdir("Logs/Damped oscillation mpi 1");
    return _execution_time_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao,
                                                      Vo, FI,
                                                      time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 * This function used to call _execution_time_damped_os_parallel_mpi_omp
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
double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
    _mkdir("Logs");
    _mkdir("Logs/Damped oscillation mpi 2");
    return _execution_time_damped_os_parallel_mpi(max_amplitude, length, mass, gravity, k, Ao,
                                                  Vo, FI,
                                                  time_limit, step_size, damping_coefficent, number_of_files);
}
/**
 * This function used to call _execution_time_damped_os_parallel_omp
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
 * @param num_of_threads number of threads needed to execute the code
 * @return
 */
double
damped_os_parallel_execution_time_v3(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files, int num_of_threads)
{
    _mkdir("Logs");
    _mkdir("Logs/Damped oscillation omp");

    return _execution_time_damped_os_parallel_omp(max_amplitude, length, mass, gravity, k,
                                                  Ao,
                                                  Vo, FI,
                                                  time_limit, step_size, damping_coefficent,
                                                  number_of_files, num_of_threads);
}
/**
 * This function used to call _execution_time_damped_os_serial
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
char *
damped_os_serial_execution(double max_amplitude, double length, double mass, double gravity, double k,
                           double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent,
                           int number_of_files)
{
    _mkdir("Logs");
    _mkdir("Logs/amped oscillation serial");

    return _execution_time_damped_os_serial(max_amplitude, length, mass, gravity, k,
                                            Ao,
                                            Vo, FI,
                                            time_limit, step_size, damping_coefficent,
                                            number_of_files);
}
/**
 * This function used to call _execution_time_elastic_pendulum
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
char *
elastic_pendulum_execution(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    _mkdir("Logs");
    _mkdir("Logs/Elastic pendulum");
    return _execution_time_elastic_pendulum(r, length, mass, gravity, k, Ao, Xo,
                                            Yo,
                                            Vo,
                                            time_limit, step_size, damping_coefficent, number_of_files);
}

////////////////heat equation 1D parallel
int heat_equation_1D_P1_MPI(double time_step, double time_limit,
                            double space_step,
                            long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 1D MPI");
    return _simulate_heat_transfer_1D_MPI(time_step, time_limit,
                                          space_step,
                                          precision);
}

int heat_equation_1D_P1_OPENMP(double time_step, double time_limit,
                               double space_step,
                               long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 1D OpenMP");
    return _simulate_heat_transfer_1D_OPENMP(time_step, time_limit,
                                             space_step,
                                             precision);
}

////////////////heat equation 2D parallel
int heat_equation_2D_P1_MPI(double time_step, double time_limit,
                            double spaceX_step, double spaceY_step,
                            long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 2D MPI");
    return _simulate_heat_transfer_2D_MPI(time_step, time_limit,
                                          spaceX_step,
                                          spaceY_step,
                                          precision);
}

int heat_equation_2D_P1_OPENMP(double time_step, double time_limit,
                               double spaceX_step, double spaceY_step,
                               long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 2D OpenMP");
    return _simulate_heat_transfer_2D_OPENMP(time_step, time_limit,
                                             spaceX_step,
                                             spaceY_step,
                                             precision);
}

////////////////heat equation 1D serial
char* heat_equation_1D_serial(double time_step, double time_limit,
                            double space_step,
                            long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 1D Serial");

    return _simulate_heat_transfer_1D_serial(time_step, time_limit,
                                             space_step,
                                             precision);
}

////////////////heat equation 2D serial
char* heat_equation_2D_serial(double time_step, double time_limit,
                            double spaceX_step, double spaceY_step,
                            long long precision)
{
    _mkdir("Simulation");
    _mkdir("Simulation/Thermo 2D Serial");
    return _simulate_heat_transfer_2D_serial(time_step, time_limit,
                                             spaceX_step,
                                             spaceY_step,
                                             precision);
}

////////////////heat equation 1D parallel execution time without i/o
double heat_equation_execution_time_1D_P1_MPI(double time_step, double time_limit,
                                              double space_step,
                                              long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/Thermo Simulation MPI 1D");
    return _execution_time_heat_transfer_1D_MPI(time_step, time_limit,
                                                space_step,
                                                precision);
}
double heat_equation_execution_time_1D_P1_OPENMP(double time_step, double time_limit,
                                                 double space_step,
                                                 long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/Thermo Simulation execution openmp 1D");
    return _execution_time_heat_transfer_1D_OPENMP(time_step, time_limit,
                                                   space_step,
                                                   precision);
}

////////////////heat equation 2D parallel execution time without i/o
double heat_equation_execution_time_2D_P1_MPI(double time_step, double time_limit,
                                              double spaceX_step, double spaceY_step,
                                              long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/Thermo Simulation mpi 2D");
    return _execution_time_heat_transfer_2D_MPI(time_step, time_limit,
                                                spaceX_step,
                                                spaceY_step,
                                                precision);
}

double heat_equation_execution_time_2D_P1_OPENMP(double time_step, double time_limit,
                                                 double spaceX_step, double spaceY_step,
                                                 long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/Thermo Simulation execution openmp 2D");
    return _execution_time_heat_transfer_2D_OPENMP(time_step, time_limit,
                                                   spaceX_step,
                                                   spaceY_step,
                                                   precision);
}

////////////////heat equation 1D serial execution time without i/o
double heat_equation_execution_time_1D_serial(double time_step, double time_limit,
                                              double space_step,
                                              long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/DThermo Simulation execution Serial 1D");
    return _execution_time_heat_transfer_1D_serial(time_step, time_limit,
                                                   space_step,
                                                   precision);
}

////////////////heat equation 2D serial execution time without i/o
double heat_equation_execution_time_2D_serial(double time_step, double time_limit,
                                              double spaceX_step, double spaceY_step,
                                              long long precision)
{
    _mkdir("Logs");
    _mkdir("Logs/Thermo Simulation execution Serial 2D");
    return _execution_time_heat_transfer_2D_serial(time_step, time_limit,
                                                   spaceX_step,
                                                   spaceY_step,
                                                   precision);
}
/**
 * This function must be called after using any collection of MPI based functions as it finalize the MPI and de-allocate the resources
 */
void finalize()
{
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Finalize();
    }
}