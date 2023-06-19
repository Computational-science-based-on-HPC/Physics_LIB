//
// Created by jghal on 6/19/2023.
//

#ifndef PHYSICS_THERMOPARA_H
#define PHYSICS_THERMOPARA_H
#include "thermoutils.h"

extern double
_get_value(struct TimeParam* time, struct SpaceParam* space, double x, double t, int precision);

extern int
_simulate_heat_transfer_1D_MPI(struct TimeParam time_param, struct SpaceParam space_param, int precision);

extern int
_simulate_heat_transfer_1D_OPENMP(struct TimeParam time_param, struct SpaceParam space_param, int precision);

extern int
_simulate_heat_transfer_1D_OPENMP_V2(struct TimeParam time_param, struct SpaceParam space_param, int precision);

extern int
_simulate_heat_transfer_1D_MPI_OPENMP(struct TimeParam time_param, struct SpaceParam space_param, int precision);

extern int
_simulate_heat_transfer_2D_MPI(struct TimeParam time_param, 
                    struct SpaceParam space_param,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
_simulate_heat_transfer_2D_OPENMP(struct TimeParam time_param, 
                    struct SpaceParam space_param,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
_simulate_heat_transfer_2D_OPENMP_V2(struct TimeParam time_param, 
                    struct SpaceParam space_param,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
_simulate_heat_transfer_2D_MPI_OPENMP(struct TimeParam time_param, 
                    struct SpaceParam space_param,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);


#endif //PHYSICS_THERMOPARA_H
