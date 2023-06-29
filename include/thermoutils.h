#ifndef PHYSICS_THERMOUTILS_H
#define PHYSICS_THERMOUTILS_H
#define ll long long
#include <stdio.h>

extern ll
_cal_num_time(double time_step, double time_limit);

extern ll
_cal_num_space(double length, double space_step);

int printmem(FILE* fileName);

void printmemsize(FILE* fileName, char *str, unsigned long ramsize);

int cpu_inf(FILE* fileName);

#endif //PHYSICS_THERMOUTILS_H