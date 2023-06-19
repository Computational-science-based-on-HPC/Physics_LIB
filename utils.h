//
// Created by jghal on 6/16/2023.
//

#ifndef PHYSICS_UTILS_H
#define PHYSICS_UTILS_H

int
_valid_osc(double x, double y, double length, double mass, double gravity, double k, double time_limit,
           double step_size,
           double damping_coefficent, int number_of_files, double Fo);

int
_min_int(int x, int y);

int
_round(double x);

double
_dx(double dx);

double
_dy(double dy);

double
_f1(double x, double y, double dx, double tx, double ty, double k, double m, double b, double r);

double
_f2(double x, double y, double dy, double tx, double ty, double k, double m, double b, double r, double g);


#endif //PHYSICS_UTILS_H
