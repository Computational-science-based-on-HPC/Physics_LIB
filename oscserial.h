//
// Created by jghal on 6/16/2023.
//

#ifndef PHYSICS_OSCSERIAL_H
#define PHYSICS_OSCSERIAL_H

extern int
_simulate_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
_simulate_forced_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                  double Vo,
                                  double time_limit, double step_size, double damping_coefficent, int number_of_files);
extern int
_simulate_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);

#endif //PHYSICS_OSCSERIAL_H
