#include <../include/physics.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
int main(void)

{
        printf("%s",damped_os_parallel_v2(10.0, 14.0, 1.0, 1.8, 1.0, -1.0, 0.0,
                         -0.1,
                         200.0, 0.1, 0.1, 3));
        finalize();
}
