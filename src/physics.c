#include "../include/physics.h"
#include "../include/oscserial.h"
#include "../include/oscpara.h"
#include "../include/thermopara.h"
#include "../include/thermoserial.h"

int damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                     double FI,
                     double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                      time_limit, step_size, damping_coefficent, number_of_files);
}

int damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                          double FI,
                          double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return _simulate_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                                time_limit, step_size, damping_coefficent, number_of_files);
}

int damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                          double FI,
                          double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return _simulate_damped_os_parallel_mpi(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                            time_limit, step_size, damping_coefficent, number_of_files);
}

int elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                     double Yo,
                     double Vo,
                     double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return _simulate_elastic_pendulum(r, length, mass, gravity, k, Ao, Xo,
                                      Yo,
                                      Vo,
                                      time_limit, step_size, damping_coefficent, number_of_files);
}

double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
    return _execution_time_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao,
                                                      Vo, FI,
                                                      time_limit, step_size, damping_coefficent, number_of_files);
}

double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
    return _execution_time_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao,
                                                      Vo, FI,
                                                      time_limit, step_size, damping_coefficent, number_of_files);
}

double
damped_os_serial_execution(double max_amplitude, double length, double mass, double gravity, double k,
                           double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent,
                           int number_of_files)
{
    return _execution_time_damped_os_serial(max_amplitude, length, mass, gravity, k,
                                            Ao,
                                            Vo, FI,
                                            time_limit, step_size, damping_coefficent,
                                            number_of_files);
}

double
elastic_pendulum_execution(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return _execution_time_elastic_pendulum(r, length, mass, gravity, k, Ao, Xo,
                                            Yo,
                                            Vo,
                                            time_limit, step_size, damping_coefficent, number_of_files);
}

////////////////heat equation 1D parallel
int heat_equation_1D_P1_MPI(double time_step, double time_limit,
                            double length, double space_step,
                            int precision)
{

    return _simulate_heat_transfer_1D_MPI(time_step, time_limit,
                                          length, space_step,
                                          precision);
}

int heat_equation_1D_P1_OPENMP(double time_step, double time_limit,
                               double length, double space_step,
                               int precision)
{

    return _simulate_heat_transfer_1D_OPENMP(time_step, time_limit,
                                             length, space_step,
                                             precision);
}

int heat_equation_1D_P1_OPENMP_V2(double time_step, double time_limit,
                                  double length, double space_step,
                                  int precision)
{

    return _simulate_heat_transfer_1D_OPENMP_V2(time_step, time_limit,
                                                length, space_step,
                                                precision);
}

////////////////heat equation 2D parallel
int heat_equation_2D_P1_MPI(double time_step, double time_limit,
                            double length, double spaceX_step, double width, double spaceY_step,
                            int precision)
{

    return _simulate_heat_transfer_2D_MPI(time_step, time_limit,
                                          length, spaceX_step,
                                          width, spaceY_step,
                                          precision);
}

int heat_equation_2D_P1_OPENMP(double time_step, double time_limit,
                               double length, double spaceX_step, double width, double spaceY_step,
                               int precision)
{

    return _simulate_heat_transfer_2D_OPENMP(time_step, time_limit,
                                             length, spaceX_step,
                                             width, spaceY_step,
                                             precision);
}

int heat_equation_2D_P1_OPENMP_V2(double time_step, double time_limit,
                                  double length, double spaceX_step, double width, double spaceY_step,
                                  int precision)
{

    return _simulate_heat_transfer_2D_OPENMP_V2(time_step, time_limit,
                                                length, spaceX_step,
                                                width, spaceY_step,
                                                precision);
}

////////////////heat equation 1D serial
int heat_equation_1D_serial(double time_step, double time_limit, double length, double space_step, int precision)
{

    return _simulate_heat_transfer_1D_serial(time_step, time_limit, length, space_step, precision);
}

////////////////heat equation 2D serial
int heat_equation_2D_serial(double time_step, double time_limit,
                            double length, double spaceX_step, double width, double spaceY_step,
                            int precision)
{

    return _simulate_heat_transfer_2D_serial(time_step, time_limit,
                                             length, spaceX_step,
                                             width, spaceY_step,
                                             precision);
}

////////////////heat equation 1D parallel execution time without i/o
int heat_equation_execution_time_1D_P1_MPI(double time_step, double time_limit,
                                           double length, double space_step,
                                           int precision)
{
    return _execution_time_heat_transfer_1D_MPI(time_step, time_limit,
                                                length, space_step,
                                                precision);
}
int heat_equation_execution_time_1D_P1_OPENMP(double time_step, double time_limit,
                                              double length, double space_step,
                                              int precision)
{
    return _execution_time_heat_transfer_1D_OPENMP(time_step, time_limit,
                                                   length, space_step,
                                                   precision);
}
int heat_equation_execution_time_1D_P1_OPENMP_V2(double time_step, double time_limit,
                                                 double length, double space_step,
                                                 int precision)
{
    return _execution_time_heat_transfer_1D_OPENMP_V2(time_step, time_limit,
                                                      length, space_step,
                                                      precision);
}

////////////////heat equation 2D parallel execution time without i/o
int heat_equation_execution_time_2D_P1_MPI(double time_step, double time_limit,
                                           double length, double spaceX_step, double width, double spaceY_step,
                                           int precision)
{
    return _execution_time_heat_transfer_2D_MPI(time_step, time_limit,
                                                length, spaceX_step,
                                                width, spaceY_step,
                                                precision);
}

int heat_equation_execution_time_2D_P1_OPENMP(double time_step, double time_limit,
                                          double length, double spaceX_step, double width, double spaceY_step,
                                          int precision){

    return _execution_time_heat_transfer_2D_OPENMP(time_step, time_limit,
                                             length, spaceX_step,
                                             width, spaceY_step,
                                             precision);
}

int heat_equation_execution_time_2D_P1_OPENMP_V2(double time_step, double time_limit,
                                             double length, double spaceX_step, double width, double spaceY_step,
                                             int precision){

    return _execution_time_heat_transfer_2D_V2_OPENMP(time_step, time_limit,
                                                length, spaceX_step,
                                                width, spaceY_step,
                                                precision);
}

////////////////heat equation 1D serial execution time without i/o
int heat_equation_execution_time_1D_serial(double time_step, double time_limit, double length, double space_step, int precision)
{
    return _execution_time_heat_transfer_1D_serial(time_step, time_limit, length, space_step, precision);
}

////////////////heat equation 2D serial execution time without i/o
int heat_equation_execution_time_2D_serial(double time_step, double time_limit,
                                           double length, double spaceX_step, double width, double spaceY_step,
                                           int precision)
{

    return _execution_time_heat_transfer_2D_serial(time_step, time_limit,
                                                   length, spaceX_step,
                                                   width, spaceY_step,
                                                   precision);
}
