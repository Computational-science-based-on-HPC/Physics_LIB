#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H

extern int
damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                 double FI,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
damped_os_parallel_v1(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
damped_os_parallel_v2(double max_amplitude, double length, double mass, double gravity, double k, double Ao, double Vo,
                      double FI,
                      double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern int
elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                 double Yo,
                 double Vo,
                 double time_limit, double step_size, double damping_coefficent, int number_of_files);

extern double
damped_os_parallel_execution_time_v1(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);
extern double
damped_os_parallel_execution_time_v2(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files);extern double
damped_os_parallel_execution_time_v3(double max_amplitude, double length, double mass, double gravity, double k,
                                     double Ao,
                                     double Vo, double FI,
                                     double time_limit, double step_size, double damping_coefficent,
                                     int number_of_files,int num_of_threads);
extern double
damped_os_serial_execution(double max_amplitude, double length, double mass, double gravity, double k,
                           double Ao,
                           double Vo, double FI,
                           double time_limit, double step_size, double damping_coefficent,
                           int number_of_files);
extern double
elastic_pendulum_execution(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                           double Yo,
                           double Vo,
                           double time_limit, double step_size, double damping_coefficent, int number_of_files);


////////////////heat equation 1D parallel
extern int
heat_equation_1D_P1_MPI(double time_step, double time_limit,
                        double space_step,
                        int precision);

extern int
heat_equation_1D_P1_OPENMP(double time_step, double time_limit,
                           double space_step,
                           int precision);

extern int
heat_equation_1D_P1_OPENMP_V2(double time_step, double time_limit,
                               double space_step,
                              int precision);

////////////////heat equation 2D parallel
extern int
heat_equation_2D_P1_MPI(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);

extern int
heat_equation_2D_P1_OPENMP(double time_step, double time_limit,
                           double spaceX_step, double spaceY_step,
                           int precision);

extern int
heat_equation_2D_P1_OPENMP_V2(double time_step, double time_limit,
                              double spaceX_step, double spaceY_step,
                              int precision);

////////////////heat equation 1D serial
extern int
heat_equation_1D_serial(double time_step, double time_limit, double space_step, int precision);


////////////////heat equation 2D serial
extern int
heat_equation_2D_serial(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);


////////////////heat equation 1D parallel execution time without i/o
extern double
heat_equation_execution_time_1D_P1_MPI(double time_step, double time_limit,
                                       double space_step,
                                       int precision);

extern double
heat_equation_execution_time_1D_P1_OPENMP(double time_step, double time_limit,
                                          double space_step,
                                          int precision);

extern double
heat_equation_execution_time_1D_P1_OPENMP_V2(double time_step, double time_limit,
                                             double space_step,
                                             int precision);

////////////////heat equation 2D parallel execution time without i/o
extern double
heat_equation_execution_time_2D_P1_MPI(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       int precision);

extern double heat_equation_execution_time_2D_P1_OPENMP(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);

extern double heat_equation_execution_time_2D_P1_OPENMP_V2(double time_step, double time_limit,
                        double spaceX_step, double spaceY_step,
                        int precision);


////////////////heat equation 1D serial execution time without i/o
extern double
heat_equation_execution_time_1D_serial(double time_step, double time_limit, double space_step, int precision);

////////////////heat equation 2D serial execution time without i/o
extern double
heat_equation_execution_time_2D_serial(double time_step, double time_limit,
                                       double spaceX_step, double spaceY_step,
                                       int precision);
extern void 
finalize();

#endif // PHYSICS_PHYSICS_H