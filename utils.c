//
// Created by jghal on 6/16/2023.
//

#include "utils.h"

int
_valid_osc(double max_amplitude, double length, double mass, double gravity, double k, double time_limit,
           double step_size,
           double damping_coefficent, int number_of_files, double Fo)
{
    if(mass<0 || k<0 || time_limit<=0 || step_size<0 || length<=0 || gravity<=0 || damping_coefficent<0 ||
       number_of_files<1)
    {
        return 0;
    }
    if(max_amplitude>length)
    {
        return -1;
    }
    return number_of_files;
}

int
_min_int(int x, int y)
{
    return x>y ? y : x;
}

int
_round(double x)
{
    return (int) (x + 0.5);
}