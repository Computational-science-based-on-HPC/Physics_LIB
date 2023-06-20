//
// Created by jghal on 6/19/2023.
//

#ifndef PHYSICS_THERMOPARA_H
#define PHYSICS_THERMOPARA_H
#include "thermoutils.h"

extern double
_get_value_1D_mpi(double time_step, 
                    double space_step, 
                    double x, double t, 
                    int precision);

// extern double
// _get_value_1D_openmp_V1(struct TimeParam time, 
//                         struct SpaceParam space, 
//                         double x, double t, 
//                         int precision);

// extern double
// _get_value_1D_openmp_V2(struct TimeParam time, 
//                         struct SpaceParam space, 
//                         double x, double t, 
//                         int precision);

// extern double
// _get_value_2D_mpi(struct TimeParam time, 
//                         struct SpaceParam2D space, 
//                         int x, int y, int t, 
//                         int precision);

// extern double
// _get_value_2D_openmp(struct TimeParam time, 
//                         struct SpaceParam2D space, 
//                         int x, int y, int t, 
//                         int precision);

// extern double
// _get_value_2D_openmp_v2(struct TimeParam time, 
//                         struct SpaceParam2D space, 
//                         int x, int y, int t, 
//                         int precision);

extern int
_simulate_heat_transfer_1D_MPI(double time_step, double time_limit, 
                        double length, double diffusivity, 
                        double space_step, 
                        int precision);

// extern int
// _simulate_heat_transfer_1D_OPENMP(struct TimeParam time_param, 
//                         struct SpaceParam space_param, 
//                         int precision);

// extern int
// _simulate_heat_transfer_1D_OPENMP_V2(struct TimeParam time_param, 
//                         struct SpaceParam space_param, 
//                         int precision);

// extern int
// _simulate_heat_transfer_1D_MPI_OPENMP(struct TimeParam time_param, 
//                         struct SpaceParam space_param, 
//                         int precision);

// extern int
// _simulate_heat_transfer_2D_MPI(struct TimeParam time_param, 
//                     struct SpaceParam2D space_param,
//                     struct TempParam temp_param, 
//                     int precision);

// extern int
// _simulate_heat_transfer_2D_OPENMP(struct TimeParam time_param, 
//                     struct SpaceParam2D space_param,
//                     struct TempParam temp_param, 
//                     int precision);

// extern int
// _simulate_heat_transfer_2D_OPENMP_V2(struct TimeParam time_param, 
//                     struct SpaceParam2D space_param,
//                     struct TempParam temp_param, 
//                     int precision);

// extern int
// _simulate_heat_transfer_2D_MPI_OPENMP(struct TimeParam time_param, 
//                     struct SpaceParam2D space_param,
//                     struct TempParam temp_param, 
//                     int precision);


#endif //PHYSICS_THERMOPARA_H
