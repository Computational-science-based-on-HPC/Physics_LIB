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
#ifdef __linux__
    return _simulate_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                                time_limit, step_size, damping_coefficent, number_of_files);
#else
    return
            _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                       time_limit, step_size, damping_coefficent, number_of_files);
#endif
}

int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
#ifdef __linux__
    return _simulate_damped_os_parallel_mpi(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                            time_limit, step_size, damping_coefficent, number_of_files);
#else
    return
            _simulate_damped_os_serial(max_amplitude, length, mass, gravity, k, Ao, Vo, FI,
                                       time_limit, step_size, damping_coefficent, number_of_files);
#endif
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


double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
#ifdef __linux__
    return _execution_time_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao,
                                                      Vo, FI,
                                                      time_limit, step_size, damping_coefficent, number_of_files);
#else
    return _execution_time_damped_os_serial(max_amplitude, length, mass, gravity, k,
                                            Ao,
                                            Vo, FI,
                                            time_limit, step_size, damping_coefficent,
                                            number_of_files);
#endif
}

double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files)
{
#ifdef __linux__
    return _execution_time_damped_os_parallel_mpi_omp(max_amplitude, length, mass, gravity, k, Ao,
                                                      Vo, FI,
                                                      time_limit, step_size, damping_coefficent, number_of_files);
#else
    return _execution_time_damped_os_serial(max_amplitude, length, mass, gravity, k,
                                            Ao,
                                            Vo, FI,
                                            time_limit, step_size, damping_coefficent,
                                            number_of_files);
#endif
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

int 
heat_equation_1D_P1_MPI(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision){
    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam space_param = {length, diffusivity, space_step};
    
    return _simulate_heat_transfer_1D_MPI(time_param, space_param, precision);
}

int
heat_equation_1D_P1_OPENMP(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision){
    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam space_param = {length, diffusivity, space_step};
    
    return _simulate_heat_transfer_1D_OPENMP(time_param, space_param, precision);
}

int
heat_equation_1D_P1_OPENMP_V2(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision){
    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam space_param = {length, diffusivity, space_step};
    
    return _simulate_heat_transfer_1D_OPENMP_V2(time_param, space_param, precision);
}

int
heat_equation_1D_P1__MPI_OPENMP(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision){
    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam space_param = {length, diffusivity, space_step};
    
    return _simulate_heat_transfer_1D_MPI_OPENMP(time_param, space_param, precision);
}

int
heat_equation_2D_P1_MPI(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision){

    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam2D space_param = {length, diffusivity, spaceX_step, width, spaceY_step};
    struct TempParam temp_param = {tempUp, tempDown, tempLeft, tempRight};
    
    return _simulate_heat_transfer_2D_MPI(time_param, space_param, temp_param, precision);
}

int
heat_equation_2D_P1_OPENMP(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision){

    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam2D space_param = {length, diffusivity, spaceX_step, width, spaceY_step};
    struct TempParam temp_param = {tempUp, tempDown, tempLeft, tempRight};
    
    return _simulate_heat_transfer_2D_OPENMP(time_param, space_param, temp_param, precision);
}

int
heat_equation_2D_P1_OPENMP_V2(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision){

    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam2D space_param = {length, diffusivity, spaceX_step, width, spaceY_step};
    struct TempParam temp_param = {tempUp, tempDown, tempLeft, tempRight};
    
    return _simulate_heat_transfer_2D_OPENMP_V2(time_param, space_param, temp_param, precision);
}

int
heat_equation_2D_P1_MPI_OPENMP(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision){

    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam2D space_param = {length, diffusivity, spaceX_step, width, spaceY_step};
    struct TempParam temp_param = {tempUp, tempDown, tempLeft, tempRight};
    
    return _simulate_heat_transfer_2D_MPI_OPENMP(time_param, space_param, temp_param, precision);
}

int heat_equation_1D_serial(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision){
    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam space_param = {length, diffusivity, space_step};
    
    return _simulate_heat_transfer_1D_serial(time_param, space_param, precision);
}

int heat_equation_2D_serial(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision){

    struct TimeParam time_param = {time_step, time_limit};
    struct SpaceParam2D space_param = {length, diffusivity, spaceX_step, width, spaceY_step};
    struct TempParam temp_param = {tempUp, tempDown, tempLeft, tempRight};
    
    return _simulate_heat_transfer_2D_serial(time_param, space_param, temp_param, precision);
}
