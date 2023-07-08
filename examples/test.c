#include <../include/physics.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
int main(void)

{

        long long precision = 2;
        for (int i = 1; i <= 10; i++)
        {
                heat_equation_execution_time_2D_serial(0.01, 0.5,0.5, 0.05, precision);

                // heat_equation_execution_time_1D_P1_MPI(0.1, 5, 0.1, 0.1, precision);
                precision *= 2;
                sleep(1);
        }
         precision = 4096;

        for (int i = 1; i <= 3; i++)
        {
                heat_equation_execution_time_2D_P1_OPENMP(0.01, 0.5,0.5, 0.05, precision);

                // heat_equation_execution_time_1D_P1_MPI(0.1, 5, 0.1, 0.1, precision);
                precision *= 2;
                sleep(1);
        }
        precision = 4096/2;
        for (int i = 1; i <= 4; i++)
        {
                heat_equation_execution_time_2D_serial(0.01, 0.5,0.5, 0.05, precision);

                // heat_equation_execution_time_1D_P1_MPI(0.1, 5, 0.1, 0.1, precision);
                precision *= 2;
                sleep(1);
        }
        finalize();
}
