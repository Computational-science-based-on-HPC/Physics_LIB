#include "physics.h"

int main(){
    // damped_os_serial(10, 1, 50, 3.8, 0.5, 0.5, 0, 0, 10, 0.01, 0.1, 10);
    heat_equation_1D_P1_MPI(0.01, 0.5, 10, 0.05, 50);
}