#include "physics.h"

int main(){
    // damped_os_serial(10, 1, 50, 3.8, 0.5, 0.5, 0, 0, 10, 0.01, 0.1, 10);
    _simulate_heat_transfer_1D_MPI(0.01, 0.5, 10, 0.1, 0.05, 50);
}