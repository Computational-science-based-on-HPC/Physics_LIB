//
// Created by jghal on 6/16/2023.
//
/**
 * @file oscserial.c
 * @brief This file contains the implementation of the serial version of the oscillation simulation in 1D and 2D.
 *
 */
#include <stdio.h>
#include "../include/oscserial.h"
#include <math.h>
#include "../include/utils.h"
#include <time.h>
/**
 *  @brief This function simulates simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @return
 */
int _simulate_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                               double Vo, double FI,
                               double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("Invalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        max_amplitude = length;
        puts("Max Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    double Wo = sqrt(k / mass);
    double W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
    double RESULTS[3];
    RESULTS[0] = max_amplitude;
    RESULTS[1] = Vo + gravity * 0;
    RESULTS[2] = Ao + gravity * exp((-damping_coefficent / (2 * mass)) * 0);
    double CALCULATIONS[3];
    FILE *p_file;
    char _file_name[2076];
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    sprintf(_file_name, "Simulation/Damped oscillation serial/damped_os_serial_displacement_%d-%02d-%02d_%02d-%02d-%02d.sim", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    p_file = fopen(_file_name, "w");
    fprintf(p_file, "%lf\n", RESULTS[0]);
    unsigned long long _it_number = (unsigned long long)((double)((time_limit / step_size) + 0.5));
    double t = 0;
    for (unsigned long long i = 0; i < _it_number; i++)
    {
        if (isnan(RESULTS[0]) || isnan(RESULTS[1]) || isnan(RESULTS[2]))
        {

            fclose(p_file);
            puts("Simulation Got a NaN Value.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a NaN Value Occurred");
            return -1;
        }
        else if (isinf(RESULTS[0]) || isinf(RESULTS[1]) || isinf(RESULTS[2]))
        {

            fclose(p_file);
            puts("Simulation Got a INF.\n Breaking the Function...\nFiles Saved.\nSimulation Ended Cause a INF Value Occurred");
            return -2;
        }
        t += step_size;
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
        fprintf(p_file, "%lf\n", RESULTS[0]);
    }
    fclose(p_file);
    return 0;
}
/**
 * @brief This function simulates the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
 *
 * @param r rest length of spring
 * @param length max length of spring
 * @param mass mass of bob suspended in spring
 * @param gravity
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return
 */
int _simulate_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                               double Yo,
                               double Vo,
                               double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    int validation = _valid_osc(Xo, Yo, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("Invalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        Xo = sqrt(length * length - Yo * Yo);
        puts("Max Amplitude Is More Than The Spring Length, Xo Had Been Reset to be Equal to the Spring Length");
    }

    double t = 0;
    double x1 = Xo; // init position of mass in x
    double y = Yo;  // init position of mass in y
    double v = Vo;  // init velocity
    double a = 0;   // init velocity
    FILE *p_dis_x, *p_dis_y;
    char _file_name[2076];
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    sprintf(_file_name, "Simulation/Elastic pendulum/elastic_pendulum_x_%d-%02d-%02d_%02d-%02d-%02d.sim", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    p_dis_x = fopen(_file_name, "w");
    sprintf(_file_name, "Simulation/Elastic pendulum/elastic_pendulum_y_%d-%02d-%02d_%02d-%02d-%02d.sim", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    p_dis_y = fopen(_file_name, "w");
    fprintf(p_dis_x, "%.6f\n", x1);
    fprintf(p_dis_y, "%.6f\n", y);

    unsigned long long _it_number = (unsigned long long)((double)((time_limit / step_size) + 0.5));
    for (unsigned long long i = 0; i < _it_number; ++i)
    {
        t += step_size;
        double k11 = _dx(v);
        double k12 = _dy(a);
        double k13 = _f1(x1, y, v, 0, 0, k, mass, damping_coefficent, r);
        double k14 = _f2(x1, y, a, 0, 0, k, mass, damping_coefficent, r, gravity);

        double k21 = _dx(v + step_size / 2 * k13);
        double k22 = _dy(a + step_size / 2 * k14);
        double k23 = _f1(x1 + step_size / 2 * k11, y + step_size / 2 * k12, v + step_size / 2 * k13, 0, 0, k, mass,
                         damping_coefficent, r);
        double k24 = _f2(x1 + step_size / 2 * k11, y + step_size / 2 * k12, a + step_size / 2 * k14, 0, 0, k, mass,
                         damping_coefficent, r, gravity);

        double k31 = _dx(v + step_size / 2 * k23);
        double k32 = _dy(a + step_size / 2 * k24);
        double k33 = _f1(x1 + step_size / 2 * k21, y + step_size / 2 * k22, v + step_size / 2 * k23, 0, 0, k, mass,
                         damping_coefficent, r);
        double k34 = _f2(x1 + step_size / 2 * k21, y + step_size / 2 * k22, a + step_size / 2 * k24, 0, 0, k, mass,
                         damping_coefficent, r, gravity);

        double k41 = _dy(v + step_size / 2 * k33);
        double k42 = _dy(a + step_size / 2 * k34);
        double k43 = _f1(x1 + step_size * k31, y + step_size * k32, v + step_size * k33, 0, 0, k, mass,
                         damping_coefficent, r);
        double k44 = _f2(x1 + step_size * k31, y + step_size * k32, a + step_size * k34, 0, 0, k, mass,
                         damping_coefficent, r, gravity);
        x1 += step_size / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        y += step_size / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        v += step_size / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        a += step_size / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);
        fprintf(p_dis_x, "%.6f\n", x1);
        fprintf(p_dis_y, "%.6f\n", y);
    }
    fclose(p_dis_x);
    fclose(p_dis_y);
    return 0;
}
/**
 *  @brief This function calculates execution time of simulating simple harmonic motion (Simple Spring Motion) using numerical solution of stepwise precision using equation (e^(-damping_coefficent / (2 * mass)) * t)*sin(wt+fi)),
 *  where this equation calculates the displacement of mass on y-axis, this function also calculates the acceleration and velocity in each time step.
 *  This function is implemented using serial algorithm
 *
 * @param max_amplitude starting position of the mass where the simulation will start
 * @param length the maximum length of the spring (uncompressed spring)
 * @param mass mass of bob
 * @param gravity
 * @param k stiffness of the spring
 * @param Ao initial acceleration
 * @param Vo initial velocity
 * @param FI FI constant which will be added to the (wt) inside the sine calculation
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files currently nulled
 * @return
 */
double
_execution_time_damped_os_serial(double max_amplitude, double length, double mass, double gravity, double k, double Ao,
                                 double Vo, double FI,
                                 double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    FILE *file;
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    char _file_name[255];
    sprintf(_file_name, "Logs/Damped oscillation serial/DOS_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    file = fopen(_file_name, "w"); // Open the file in write mode
    if (file != NULL)
    {
        freopen(_file_name, "w", stdout);
    }
    printf("Started Simulation of Damped Oscillation Serial Implementation at %d-%02d-%02d %02d:%02d:%02d with Parameters:\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    printf("Amplitude: %f \nSpring Length: %f \nMass: %f \nGravity: %f \nStifeness: %f \nInitial Acceleration: %f \nInitial Velocity: %f \nFI Const: %f \nTime Limit: %f \nStep_Size(dt): %f \nDamping coefficient: %f\n", max_amplitude, length, mass, gravity, k, Ao, Vo, FI, time_limit, step_size, damping_coefficent);
    puts("\n================================================================================");
    printf("Memory ===========================================================================\n");
    puts("================================================================================");
    printmemstream();
    printf("\n================================================================================\nCPUs ===========================================================================\n");
    puts("================================================================================\n\n");
    cpu_inf_stream();
    printf("\n================================================================================\nExecution Times===========================================================================\n");
    puts("================================================================================\n\n");
    int validation = _valid_osc(max_amplitude, 0, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("\n\nInvalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        max_amplitude = length;
        puts("\n\nMax Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    double Wo = sqrt(k / mass);
    double W = sqrt(Wo - pow(damping_coefficent / 2 * mass, 2));
    double RESULTS[3];
    RESULTS[0] = max_amplitude;
    RESULTS[1] = Vo + gravity * 0;
    RESULTS[2] = Ao + gravity * exp((-damping_coefficent / (2 * mass)) * 0);
    double CALCULATIONS[3];
    tm = *localtime(&tim);
    printf("\nStarted Calculation at: %d-%02d-%02d %02d:%02d:%02d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    clock_t start_time = clock();
    unsigned long long _it_number = (unsigned long long)((double)((time_limit / step_size) + 0.5));
    double t = 0;
    for (unsigned long long i = 0; i < _it_number; i++)
    {
        t += step_size;
        CALCULATIONS[0] = cos(W * t + FI);
        CALCULATIONS[1] = sin(W * t + FI);
        CALCULATIONS[2] = exp((-damping_coefficent / (2 * mass)) * t);
        RESULTS[0] = max_amplitude * CALCULATIONS[2] * CALCULATIONS[0];
        RESULTS[1] = W * max_amplitude * CALCULATIONS[2] * CALCULATIONS[1] + (gravity * t * CALCULATIONS[2]);
        RESULTS[2] = (-1 * W * W * CALCULATIONS[2] * CALCULATIONS[0]) + (gravity * CALCULATIONS[2]);
    }
    clock_t end_time = clock();
    tim = time(NULL);
    tm = *localtime(&tim);
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\nEnded Job Caclulation at: %d-%02d-%02d %02d:%02d:%02d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, execution_time);
    puts("\n================================================================================\n");
    printf("Ended Simulation of Damped Oscillation Implementation at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    return execution_time;
}
/**
 * @brief This function calculates the execution time of simulating the motion of (elastic pendulum/2D-spring/spring pendulum) system.
 * Using LaGrange mechanics to get the equation of motion of the whole system and solving the differential equation using Fourth order Runge-Kutta ODE to get the displacement of body suspended on spring at time t.
 * This system's motion is chaotic motion so it can't be parallelized.
 * This simulation prints the position of mass w.r.t X-Axis and Y-Axis.
 *
 * @param r rest length of spring
 * @param length max length of spring
 * @param mass mass of bob suspended in spring
 * @param gravity
 * @param k stiffness of spring
 * @param Ao initial acceleration
 * @param Xo initial point on X-axis where simulation starts
 * @param Yo initial point on Y-axis where simulation starts
 * @param Vo initial velocity
 * @param time_limit the time when the simulation will stop
 * @param step_size how much the simulation will skip per iteration
 * @param damping_coefficent damping factor affecting on the system
 * @param number_of_files
 * @return
 */
double
_execution_time_elastic_pendulum(double r, double length, double mass, double gravity, double k, double Ao, double Xo,
                                 double Yo,
                                 double Vo,
                                 double time_limit, double step_size, double damping_coefficent, int number_of_files)
{
    time_t tim = time(NULL);
    struct tm tm = *localtime(&tim);
    FILE *file;
    char _file_name[255];
    sprintf(_file_name, "Logs/Elastic Pendulum/EP_%d-%02d-%02d_%02d-%02d-%02d.log", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    file = fopen(_file_name, "w"); // Open the file in write mode
    if (file != NULL)
    {
        freopen(_file_name, "w", stdout);
    }
    printf("Started Simulation of Elastic Pendulum Oscillation Serial Implementation at %d-%02d-%02d %02d:%02d:%02d with Parametes:\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    printf("Rest Length: %f \nSpring Length: %f \nMass: %f \nGravity: %f \nStifeness: %f \nInitial Acceleration: %f \nInitial Velocity: %f \nInitial Y: %f\nInitial X: %f\nTime Limit: %f \nStep_Size(dt): %f \nDamping_coefficient: %f\n", r, length, mass, gravity, k, Ao, Vo, Yo, Xo, time_limit, step_size, damping_coefficent);
    puts("\n================================================================================");
    printf("Memory ===========================================================================\n");
    puts("================================================================================");
    printmemstream();
    printf("\n================================================================================\nCPUs ===========================================================================\n");
    puts("================================================================================\n\n");
    cpu_inf_stream();
    printf("\n================================================================================\nExecution Times===========================================================================\n");
    puts("================================================================================\n\n");
    int validation = _valid_osc(Xo, Yo, length, mass, gravity, k, time_limit, step_size, damping_coefficent,
                                number_of_files, 0);
    if (validation == 0)
    {
        puts("\n\nInvalid Arguments is Given");
        return -1;
    }
    if (validation == -1)
    {
        Xo = sqrt(length * length - Yo * Yo);
        puts("\n\nMax Amplitude Is More Than The Spring Length, Max Amplitude is Set Equal to Spring Length");
    }
    double t = 0;
    double x1 = Xo; // init position of mass in x
    double y = Yo;  // init position of mass in y
    double v = Vo;  // init velocity
    double a = 0;   // init velocity
    unsigned long long _it_number = (unsigned long long)((double)((time_limit / step_size) + 0.5));
    tm = *localtime(&tim);
    printf("\nStarted Calculation at: %d-%02d-%02d %02d:%02d:%02d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    clock_t start_time = clock();
    for (unsigned long long i = 0; i < _it_number; ++i)
    {
        t += step_size;
        double k11 = _dx(v);
        double k12 = _dy(a);
        double k13 = _f1(x1, y, v, 0, 0, k, mass, damping_coefficent, r);
        double k14 = _f2(x1, y, a, 0, 0, k, mass, damping_coefficent, r, gravity);

        double k21 = _dx(v + step_size / 2 * k13);
        double k22 = _dy(a + step_size / 2 * k14);
        double k23 = _f1(x1 + step_size / 2 * k11, y + step_size / 2 * k12, v + step_size / 2 * k13, 0, 0, k, mass,
                         damping_coefficent, r);
        double k24 = _f2(x1 + step_size / 2 * k11, y + step_size / 2 * k12, a + step_size / 2 * k14, 0, 0, k, mass,
                         damping_coefficent, r, gravity);

        double k31 = _dx(v + step_size / 2 * k23);
        double k32 = _dy(a + step_size / 2 * k24);
        double k33 = _f1(x1 + step_size / 2 * k21, y + step_size / 2 * k22, v + step_size / 2 * k23, 0, 0, k, mass,
                         damping_coefficent, r);
        double k34 = _f2(x1 + step_size / 2 * k21, y + step_size / 2 * k22, a + step_size / 2 * k24, 0, 0, k, mass,
                         damping_coefficent, r, gravity);

        double k41 = _dy(v + step_size / 2 * k33);
        double k42 = _dy(a + step_size / 2 * k34);
        double k43 = _f1(x1 + step_size * k31, y + step_size * k32, v + step_size * k33, 0, 0, k, mass,
                         damping_coefficent, r);
        double k44 = _f2(x1 + step_size * k31, y + step_size * k32, a + step_size * k34, 0, 0, k, mass,
                         damping_coefficent, r, gravity);
        x1 += step_size / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        y += step_size / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        v += step_size / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        a += step_size / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);
    }
    clock_t end_time = clock();
    tim = time(NULL);
    tm = *localtime(&tim);
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\nEnded Job Caclulation at: %d-%02d-%02d %02d:%02d:%02d Execution Time: %f sec.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, execution_time);
    puts("\n================================================================================\n");
    printf("Ended Simulation of Elastic Pendulum Implementation at %d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    return execution_time;
}
