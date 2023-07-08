#include <../include/physics.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
int main(void)

{

    // long long precision = 2;
    // for (int i = 1; i <= 14; i++)
    // {
    //     heat_equation_execution_time_2D_P1_MPI(0.01, 0.5, 0.5, 0.05, precision);
    //     precision *= 2;
    //     sleep(1);
    // }
    /**
     * @brief Amplitude: 10.000000
Spring Length: 14.000000
Mass: 1.000000
Gravity: 1.800000
Stifeness: 1.000000
Initial Acceleration: -1.000000
Initial Velocity: 0.000000
FI Const: -0.100000
Time Limit: 33554432.000000
Step_Size(dt): 0.010000
Damping coefficient: 0.100000
Numper of Processes: 10
     *
     */
    // int time = 33554432;
    // for (int i = 1; i <= 3; i++)
    // {
    //     damped_os_parallel_execution_time_v2(10, 14.000000, 1, 1.8, 1, -1, 0, -0.1, time, 0.01, 0.1, 10);
    //     time *= 2;
    //     sleep(1);
    // }
    // finalize();
    int time = 1048576*2;
    for (int i = 1; i <= 10; i++)
    {
        elastic_pendulum(10, 14.000000, 1, 9.8, 1, -1, 0, -0.1, 40, 0.01, 0.1, 10,1);

        damped_os_parallel_v3(10, 14.000000,1, 1.8, 1, -1, 0, -0.1, time, 0.01, 0.1, 10,16);

        time *= 2;
        sleep(1);
    }
}
