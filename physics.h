#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H

extern int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                          double time_limit, double step_size, double damping_coefficent, int number_of_files);
extern int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
extern int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
#endif //PHYSICS_PHYSICS_H
extern int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);