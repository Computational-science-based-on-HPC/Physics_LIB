#include "physics.h"
#include "oscserial.h"
#include "oscpara.h"
#include "thermopara.h"

int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                 double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return
            _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                       time_limit, step_size, damping_coefficent, number_of_files);
}

int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
#ifdef _WIN32
    return
            _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                       time_limit, step_size, damping_coefficent, number_of_files);
#endif
    return _simulate_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                                time_limit, step_size, damping_coefficent, number_of_files);
}

int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
#ifdef _WIN32
    return
            _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                       time_limit, step_size, damping_coefficent, number_of_files);
#endif
    return _simulate_damped_os_parallel_mpi(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                            time_limit, step_size, damping_coefficent, number_of_files);
}

int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    return
            _simulate_elastic_pendulum(r, length, mass, gravity, k, Ao, Xo,
                                       Yo,
                                       Vo,
                                       time_limit, step_size, damping_coefficent, number_of_files);
}

int
heat_equation_1D_V1_MPI(double time_step, double time_limit, double length, double diffusivity, double space_step,
                        int precision)
{
    struct TimeParam time_param={time_step, time_limit};
    struct SpaceParam space_param={length, diffusivity, space_step};

    return _simulate_heat_transfer_1D_MPI(time_param, space_param, precision);
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
