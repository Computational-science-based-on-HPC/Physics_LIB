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
