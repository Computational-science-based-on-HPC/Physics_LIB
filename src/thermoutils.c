#include "thermoutils.h"
#define ll long long


void _cal_num_time(struct TimeParam* t, ll* numTimePoint){
    numTimePoint = (ll)(t->t_lim / t->delta_t);
}

void _cal_num_space(struct SpaceParam* s, ll* numSpacePoint){
    numSpacePoint = (ll)(s->length / s->delta_x);
}

void _cal_num_space_2D(struct TimeParam* time_param, struct SpaceParam2D* space_param, ll* numSpacePointX, ll* numSpacePointY){
    numSpacePointX = (ll)(space_param->length / space_param->delta_x);
    numSpacePointY = (ll)(space_param->width / space_param->delta_y);

}
