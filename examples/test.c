#include <../include/physics.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
int main(void)

{
        long long precision = 8192*2;

        for (int i = 1; i <= 10; i++)
        {
                heat_equation_execution_time_2D_P1_OPENMP(0.01, 0.5,0.5, 0.05, precision);

                // heat_equation_execution_time_1D_P1_MPI(0.1, 5, 0.1, 0.1, precision);
                precision *= 2;
                sleep(5);
        }
        finalize();
}
