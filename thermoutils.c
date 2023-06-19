#include "thermoutils.h"
#define ll long long

// struct TimeParam
// {
//     double delta_t, t_lim;
// };

// struct SpaceParam
// {
//     double length, diffusion, delta_x;
// };

void _cal_num_time(struct TimeParam* t, ll* numTimePoint){
    numTimePoint = (ll)(t->t_lim / t->delta_t);
}

void _cal_num_space(struct SpaceParam* s, ll* numSpacePoint){
    numSpacePoint = (ll)(s->length / s->delta_x);
}