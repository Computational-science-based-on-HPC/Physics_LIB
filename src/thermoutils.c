#include "../include/thermoutils.h"
#include <stdio.h>
#include <sys/sysinfo.h>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

#define ll long long

unsigned ll _cal_num_time(double time_step, double time_limit){
    return (ll)(time_limit / time_step);
}

ll _cal_num_space(double length, double space_step){
    return (ll)(length / space_step);
}
