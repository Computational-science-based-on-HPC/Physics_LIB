//
// Created by jghal on 6/16/2023.
//
#include <stdio.h>
#include "oscserial.h"
#include "math.h"
#include "utils.h"

const char *FILES[]={"displacement.txt", "velocity.txt", "acceleration.txt"};

int
_simulate_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int validation=_valid_osc(max_amplitude, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                              number_of_files, 0);
    if(validation==0)
    {
        puts("Invalid Arguments is Given");
        return -1;
    }
    if(validation==-1)
    {
        max_amplitude=length;
        puts("Max Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    double Wo=sqrt(k/mass);
    double W=sqrt(Wo - pow(damping_coefficent/2*mass, 2));
    double RESULTS[3];
    RESULTS[0]=max_amplitude;
    RESULTS[1]=Vo + gravity*0;
    RESULTS[2]=Ao + gravity*exp((-damping_coefficent/(2*mass))*0);
    double CALCULATIONS[3];
    number_of_files=_min_int(number_of_files, 3);
    FILE *p_file[number_of_files];
    for(int i=0; i<number_of_files; ++i)
    {
        p_file[i]=fopen(FILES[i], "w");
    }
    for(double t=0; t<=time_limit + 0.05; t+=step_size)
    {
        if(isnan(RESULTS[0]) || isnan(RESULTS[1]) || isnan(RESULTS[2]))
        {
            for(int i=0; i<number_of_files; ++i)
            {
                fclose(p_file[i]);
            }
            puts("Simulation Got a NaN Value.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a NaN Value Occurred");
            return -1;
        } else if(isinf(RESULTS[0]) || isinf(RESULTS[1]) || isinf(RESULTS[2]))
        {
            for(int i=0; i<number_of_files; ++i)
            {
                fclose(p_file[i]);
            }
            puts("Simulation Got a INF.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a INF Value Occurred");
            return -2;
        }
        for(int i=0; i<number_of_files; ++i)
        {
            fprintf(p_file[i], "%lf\n", RESULTS[i]);
        }
        CALCULATIONS[0]=cos(W*t + FI);
        CALCULATIONS[1]=sin(W*t + FI);
        CALCULATIONS[2]=exp((-damping_coefficent/(2*mass))*t);
        RESULTS[0]=max_amplitude*CALCULATIONS[2]*CALCULATIONS[0];
        RESULTS[1]=W*max_amplitude*CALCULATIONS[2]*CALCULATIONS[1] + (gravity*t*CALCULATIONS[2]);
        RESULTS[2]=(-1*W*W*CALCULATIONS[2]*CALCULATIONS[0]) + (gravity*CALCULATIONS[2]);

    }
    for(int i=0; i<number_of_files; ++i)
    {
        fclose(p_file[i]);
    }
    return 0;
}


int
_simulate_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    double t=0;
    double x1=Xo; // init position of mass in x
    double x2=Yo; // init position of mass in y
    double x3=Vo;   // init velocity
    double x4=0;   // init velocity
    for(int i=0; i<time_limit; ++i)
    {
        t+=step_size;
        double k11=_dx(x3);
        double k12=_dy(x4);
        double k13=_f1(x1, x2, x3, 0, 0, k, mass, damping_coefficent, r);
        double k14=_f2(x1, x2, x4, 0, 0, k, mass, damping_coefficent, r, gravity);

        double k21=_dx(x3 + step_size/2*k13);
        double k22=_dy(x4 + step_size/2*k14);
        double k23=_f1(x1 + step_size/2*k11, x2 + step_size/2*k12, x3 + step_size/2*k13, 0, 0, k, mass,
                       damping_coefficent, r);
        double k24=_f2(x1 + step_size/2*k11, x2 + step_size/2*k12, x4 + step_size/2*k14, 0, 0, k, mass,
                       damping_coefficent, r, gravity);

        double k31=_dx(x3 + step_size/2*k23);
        double k32=_dy(x4 + step_size/2*k24);
        double k33=_f1(x1 + step_size/2*k21, x2 + step_size/2*k22, x3 + step_size/2*k23, 0, 0, k, mass,
                       damping_coefficent, r);
        double k34=_f2(x1 + step_size/2*k21, x2 + step_size/2*k22, x4 + step_size/2*k24, 0, 0, k, mass,
                       damping_coefficent, r, gravity);

        double k41=_dy(x3 + step_size/2*k33);
        double k42=_dy(x4 + step_size/2*k34);
        double k43=_f1(x1 + step_size*k31, x2 + step_size*k32, x3 + step_size*k33, 0, 0, k, mass,
                       damping_coefficent, r);
        double k44=_f2(x1 + step_size*k31, x2 + step_size*k32, x4 + step_size*k34, 0, 0, k, mass,
                       damping_coefficent, r, gravity);
        x1+=step_size/6.0*(k11 + 2*k21 + 2*k31 + k41);
        x2+=step_size/6.0*(k12 + 2*k22 + 2*k32 + k42);
        x3+=step_size/6.0*(k13 + 2*k23 + 2*k33 + k43);
        x4+=step_size/6.0*(k14 + 2*k24 + 2*k34 + k44);
    }
}