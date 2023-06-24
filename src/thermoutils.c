#include "thermoutils.h"
#define ll long long


ll _cal_num_time(double time_step, double time_limit){
    // ll time = (time_limit / time_step);
    // numTimePoint = &time;
    return (ll)(time_limit / time_step);
}

ll _cal_num_space(double length, double space_step){
    // ll space = (length / space_step);
    // numSpacePoint = &space;
    return (ll)(length / space_step);
}

// void _cal_num_space_2D(struct TimeParam* time_param, struct SpaceParam2D* space_param, ll* numSpacePointX, ll* numSpacePointY){
//     numSpacePointX = (ll)(space_param->length / space_param->delta_x);
//     numSpacePointY = (ll)(space_param->width / space_param->delta_y);

// }
