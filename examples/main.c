#include "../include/physics.h"

int main(void)
{

    damped_os_serial(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
                     -0.1,
                     400, 0.01, 0.1, 3);
}