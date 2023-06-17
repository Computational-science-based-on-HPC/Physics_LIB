//
// Created by jghal on 6/16/2023.
//

#ifndef PHYSICS_UTILS_H
#define PHYSICS_UTILS_H

extern int
_valid_osc(double max_amplitude, double length, double mass, double gravity, double k, double time_limit,
           double step_size,
           double damping_coefficent, int number_of_files, double Fo);

int
_min_int(int x, int y);
int
_round(double x);


#endif //PHYSICS_UTILS_H
