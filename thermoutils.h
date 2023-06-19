

#ifndef PHYSICS_THERMOUTILS_H
#define PHYSICS_THERMOUTILS_H
#define ll long long

extern void
_cal_num_time(struct TimeParam* t, ll* numTimePoint);

extern void
_cal_num_space(struct SpaceParam* s, ll* numSpacePoint);

struct TimeParam
{
    double delta_t, t_lim;
};

struct SpaceParam
{
    double length, diffusion, delta_x;
};

#endif //PHYSICS_THERMOUTILS_H