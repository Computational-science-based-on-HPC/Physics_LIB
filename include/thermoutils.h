

#ifndef PHYSICS_THERMOUTILS_H
#define PHYSICS_THERMOUTILS_H
#define ll long long

extern void
_cal_num_time(struct TimeParam* t, ll* numTimePoint);

extern void
_cal_num_space(struct SpaceParam* s, ll* numSpacePoint);

extern void
_cal_num_space_2D(struct TimeParam* time_param, struct SpaceParam2D* space_param, ll* numSpacePointX, ll* numSpacePointY);


struct TimeParam
{
    double delta_t, t_lim;
};

struct SpaceParam
{
    double length, diffusion, delta_x;
};

struct SpaceParam2D
{
    double length, diffusion, delta_x, width, delta_y;
};

struct TempParam
{
    double tempUp, tempDown, tempLeft, tempRight;
};

#endif //PHYSICS_THERMOUTILS_H