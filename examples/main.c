#include "physics.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

int main(void)

{
//    int t = 2;
//    for (int i = 0; i < 30; i++)
//    {
//        damped_os_parallel_execution_time_v1(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
//                                             -0.1,
//                                             t, 0.01, 0.1, 3);
//        t *= 2;
//    }
//    t = 2;
//    for (int i = 0; i < 30; i++)
//    {
//        damped_os_parallel_execution_time_v2(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
//                                             -0.1,
//                                             t, 0.01, 0.1, 3);
//        t *= 2;
//    }
    //    damped_os_parallel_execution_time_v3(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
    //                                         -0.1,
    //                                         40, 0.01, 0.1, 3,10);
    // damped_os_serial_execution(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
    //                            -0.1,
    //                            40000000, 0.01, 0.1, 3);

    // damped_os_serial_execution(10.0, 14.0, 1.0, 1.8,  1.0,
    //                            -1.0,
    //                            0.0, -0.1,
    //                            40000000, 0.01, 0.1,
    //                           3);

    //    heat_equation_1D_P1_MPI(0.01, 0.5, 0.05, 50);
//        heat_equation_1D_P1_OPENMP(0.01, 0.5, 0.05, 50);
    // for(int i=1; i<= 30; i++){
    //  long long precision = 50;
    //  heat_equation_2D_P1_MPI(0.1, 5, 0.1, 0.1, precision);
    //     precision*=2;
    // }

//        heat_equation_2D_P1_OPENMP(0.1, 5, 0.1, 0.1, 50);

        heat_equation_1D_serial(0.01, 0.5, 0.05,50);
     heat_equation_2D_serial(0.1, 5, 0.1, 0.1,50);

    //    heat_equation_execution_time_1D_P1_MPI(0.01, 0.5, 0.05, 50);
//        heat_equation_execution_time_1D_P1_OPENMP(0.01, 0.5, 0.05, 50);

//    heat_equation_execution_time_2D_P1_MPI(0.1, 5, 0.1, 0.1, 50);
//    heat_equation_execution_time_2D_P1_OPENMP(0.1, 5, 0.1, 0.1, 50);

        heat_equation_execution_time_1D_serial(0.01, 0.5, 0.05,50);
     heat_equation_execution_time_2D_serial(0.1, 5, 0.1, 0.1,50);

    //  finalize();
}
