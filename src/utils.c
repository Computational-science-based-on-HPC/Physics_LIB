//
// Created by jghal on 6/16/2023.
//

#include "../include/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <sys/sysinfo.h>
#define _GNU_SOURCE
int cpu_inf_stream()
{
    FILE *cpuinfo = fopen("/proc/cpuinfo", "rb");
    char *arg = 0;
    size_t size = 0;
    while (getdelim(&arg, &size, 0, cpuinfo) != -1)
    {
        puts(arg);
    }
    free(arg);
    fclose(cpuinfo);
    return 0;
}
int _valid_osc(double x, double y, double length, double mass, double gravity, double k, double time_limit,
               double step_size,
               double damping_coefficent, int number_of_files, double Fo, double rest_length)
{
    double max_length = sqrt((x * x) + (y * y));
    if (mass < 0 || k < 0 || time_limit <= 0 || step_size < 0 || length <= 0 || gravity <= 0 || damping_coefficent < 0 || rest_length > length || Fo < 0)
    {
        return 0;
    }
    if (max_length > length)
    {
        return -1;
    }
    return 1;
}

int _min_int(int x, int y)
{
    return x > y ? y : x;
}

int _round(double x)
{
    return (int)(x + 0.5);
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
    double L = sqrt(pow(x - tx, 2) + pow(y - ty, 2));
    double s = (x - tx) / L;
    return -(double)k / m * r * s - (float)b / m * dx;
}

double
_f2(double x, double y, double dy, double tx, double ty, double k, double m, double b, double r, double g)
{
    {
        double L = sqrt(pow(x - tx, 2) + pow(y - ty, 2));
        double c = (y - ty) / L;

        return g - (float)k / m * r * c - (float)b / m * dy;
    }
}

void printmemsizestream(char *str, unsigned long ramsize)
{
    printf("%s: %ld in bytes / %ld in KB / %ld in MB / %ld in GB\n", str, ramsize, ramsize / 1024, (ramsize / 1024) / 1024, ((ramsize / 1024) / 1024) / 1024);
}

int printmemstream()
{
    struct sysinfo info;
    sysinfo(&info);
    printf("\n\nuptime: %ld\n", info.uptime);
    // print total ram size
    printmemsizestream("totalram", info.totalram);
    printmemsizestream("freeram", info.freeram);
    printmemsizestream("sharedram", info.sharedram);
    printmemsizestream("bufferram", info.bufferram);
    printmemsizestream("freeswap", info.freeswap);
    printf("current running processes: %d\n\n", info.procs);
    return 0;
}
int _mkdir(char *_dir_name)
{
    if (mkdir(_dir_name, 0777) == -1)
    {
        puts("Directory already exists");
        return -1;
    }
    else
        puts("Directory created");
    return 0;
}