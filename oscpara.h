//
// Created by jghal on 6/16/2023.
//

#ifndef PHYSICS_OSCPARA_H
#define PHYSICS_OSCPARA_H

extern int
_simulate_damped_os_parallel_mpPi_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);

extern int
_execution_time_damped_os_parallel_mpi_omp(double max_amplitude, double length, double mass, double gravity, double k,
                                           double Ao,
                                           double Vo, double FI,
                                           double time_limit, double step_size, double damping_coefficent,
                                           int number_of_files);

extern int
_simulate_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                 double Vo, double FI,
                                 double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
_execution_time_damped_os_parallel_mpi(double max_amplitude, double length, double mass, double gravity, double k,
                                       double Ao,
                                       double Vo, double FI,
                                       double time_limit, double step_size, double damping_coefficent,
                                       int number_of_files);

extern int _min_int(int x, int y);

#endif //PHYSICS_OSCPARA_H
