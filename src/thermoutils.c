#include "../include/thermoutils.h"
#include <stdio.h>
#include <sys/sysinfo.h>

#define ll long long

ll _cal_num_time(double time_step, double time_limit){
    return (ll)(time_limit / time_step);
}

ll _cal_num_space(double length, double space_step){
    return (ll)(length / space_step);
}

void printmemsize(FILE* fileName, char *str, unsigned long ramsize) {
    
    fprintf(fileName, "%s: %ld in bytes / %ld in KB / %ld in MB / %ld in GB\n",str, ramsize, ramsize/1024, (ramsize/1024)/1024, ((ramsize/1024)/1024)/1024);
}

int printmem(FILE* fileName) {
    struct sysinfo info;
    sysinfo(&info);
    fprintf(fileName, "\n\nuptime: %ld\n", info.uptime);
    // print total ram size
    printmemsize(fileName, "totalram", info.totalram);
    printmemsize(fileName, "freeram", info.freeram);
    printmemsize(fileName, "sharedram", info.sharedram);
    printmemsize(fileName, "bufferram", info.bufferram);
    printmemsize(fileName, "freeswap", info.freeswap);
    fprintf(fileName, "current running processes: %d\n\n", info.procs);
    return 0;
}

int cpu_inf(FILE* fileName)
{
   FILE *cpuinfo = fopen("/proc/cpuinfo", "rb");
   char *arg = 0;
   size_t size = 0;
   while(getdelim(&arg, &size, 0, cpuinfo) != -1)
   {
        puts(arg);
        fprintf(fileName, "%s\n", arg);
   }
   free(arg);
   fclose(cpuinfo);
   return 0;
}
