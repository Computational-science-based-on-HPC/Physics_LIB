#include "../include/thermoutils.h"
#define ll long long


ll _cal_num_time(double time_step, double time_limit){
    return (ll)(time_limit / time_step);
}

ll _cal_num_space(double length, double space_step){
    return (ll)(length / space_step);
}

