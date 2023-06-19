#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H

extern int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                          double time_limit, double step_size, double damping_coefficent, int number_of_files);
extern int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);
extern int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
heat_equation_1D_P1_MPI(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision);    

extern int
heat_equation_1D_P1_OPENMP(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision);

extern int
heat_equation_1D_P1_OPENMP_V2(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision);

extern int
heat_equation_1D_P1__MPI_OPENMP(double time_step, double time_limit, double length, double diffusivity, double space_step, int precision);

extern int
heat_equation_2D_P1_MPI(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
heat_equation_2D_P1_OPENMP(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
heat_equation_2D_P1_OPENMP_V2(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
heat_equation_2D_P1_MPI_OPENMP(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

extern int
heat_equation_1D_serial(double time_step, double time_limit, 
                    double length, double diffusivity, double space_step, 
                    int precision);

extern int
heat_equation_2D_serial(double time_step, double time_limit, 
                    double length, double diffusivity, double spaceX_step, double width, double spaceY_step,
                    double tempUp, double tempDown, double tempLeft, double tempRight, int precision);

#endif //PHYSICS_PHYSICS_H