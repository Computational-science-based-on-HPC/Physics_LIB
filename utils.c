//
// Created by jghal on 6/16/2023.
//

#include "utils.h"
#include "math.h"

int
_valid_osc(double x, double y, double length, double mass, double gravity, double k, double time_limit,
           double step_size,
           double damping_coefficent, int number_of_files, double Fo)
{
    double max_length=sqrt((x*x) + (y*y));
    if(mass<0 || k<0 || time_limit<=0 || step_size<0 || length<=0 || gravity<=0 || damping_coefficent<0 ||
       number_of_files<1)
    {
        return 0;
    }
    if(max_length>length)
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


double
_dx(double dx)
{
    return dx;
}

double
_dy(double dy)
{
    return dy;
}

double
_f1(double x, double y, double dx, double tx, double ty, double k, double m, double b, double r)
{
    double L=sqrt(pow(x - tx, 2) + pow(y - ty, 2));
    double s=(x - tx)/L;
    return -(double) k/m*r*s - (float) b/m*dx;

}

double
_f2(double x, double y, double dy, double tx, double ty, double k, double m, double b, double r, double g)
{
    {
        double L=sqrt(pow(x - tx, 2) + pow(y - ty, 2));
        double c=(y - ty)/L;

        return g - (float) k/m*r*c - (float) b/m*dy;
    }
}