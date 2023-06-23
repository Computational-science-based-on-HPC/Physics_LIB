#include "physics.h"

int main(){
    // damped_os_serial(10, 1, 50, 3.8, 0.5, 0.5, 0, 0, 10, 0.01, 0.1, 10);
//    heat_equation_1D_P1_MPI(0.01, 0.5, 10, 0.05, 50);
//    heat_equation_1D_P1_OPENMP(0.01, 0.5, 10, 0.05, 50);
//    heat_equation_1D_P1_OPENMP_V2(0.01, 0.5, 10, 0.05, 50);

   heat_equation_2D_P1_MPI(0.1, 5,2, 0.1, 2, 0.1,50);
//    heat_equation_2D_P1_OPENMP(0.1, 5,2, 0.1, 2, 0.1,50);
//    heat_equation_2D_P1_OPENMP_V2(0.1, 5,2, 0.1, 2, 0.1,50);

//    heat_equation_1D_serial(0.01, 0.5, 10, 0.05,50);
//    heat_equation_2D_serial(0.1, 5,2, 0.1, 2, 0.1,50);

//    heat_equation_execution_time_1D_P1_MPI(0.01, 0.5, 10, 0.05, 50);
    // heat_equation_execution_time_1D_P1_OPENMP(0.01, 0.5, 10, 0.05, 50);
    // heat_equation_execution_time_1D_P1_OPENMP_V2(0.01, 0.5, 10, 0.05, 50);

    heat_equation_execution_time_2D_P1_MPI(0.1, 5,2, 0.1, 2, 0.1,50);

//    heat_equation_execution_time_1D_serial(0.01, 0.5, 10, 0.05,50);
//    heat_equation_execution_time_2D_serial(0.1, 5,2, 0.1, 2, 0.1,50);


}