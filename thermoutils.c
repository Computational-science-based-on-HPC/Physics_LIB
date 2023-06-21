#include "thermoutils.h"
#define ll long long


void _cal_num_time(double time_step, double time_limit, ll* numTimePoint){
    ll time = (time_step / time_limit);
    numTimePoint = &time;
}

void _cal_num_space(double length, double space_step, ll* numSpacePoint){
    ll space = (length / space_step);
    numSpacePoint = &space;
}

// void _cal_num_space_2D(struct TimeParam* time_param, struct SpaceParam2D* space_param, ll* numSpacePointX, ll* numSpacePointY){
//     numSpacePointX = (ll)(space_param->length / space_param->delta_x);
//     numSpacePointY = (ll)(space_param->width / space_param->delta_y);

// }
