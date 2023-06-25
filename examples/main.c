#include "../include/physics.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

int main(void)
{
    damped_os_serial(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
                          -0.1,
                          40, 1, 0.1, 3);
    // finalize();
}