#ifndef PHYSICS_THERMOSERIAL_H
#define PHYSICS_THERMOSERIAL_H
#include "thermoutils.h"

extern int
_simulate_heat_transfer_1D_serial(struct TimeParam time_param, 
                        struct SpaceParam space_param, 
                        int precision);

extern double
_get_value_1D(struct TimeParam time, 
                        struct SpaceParam space, 
                        double x, double t, 
                        int precision);

extern double
_get_value_2D(struct TimeParam time, 
                        struct SpaceParam2D space, 
                        int x, int y, int t, 
                        int precision);

extern int
_simulate_heat_transfer_2D_serial(struct TimeParam time_param, 
                        struct SpaceParam2D space_param,
                        struct TempParam temp_param, 
                        int precision);

#endif //PHYSICS_THERMOSERIAL_H
