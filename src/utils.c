//
// Created by jghal on 6/16/2023.
//

#include "../include/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
int 
_valid_osc(double x, double y, double length, double mass, double gravity, double k, double time_limit,
               double step_size,
               double damping_coefficent, int number_of_files, double Fo)
{
    double max_length = sqrt((x * x) + (y * y));
    if (mass < 0 || k < 0 || time_limit <= 0 || step_size < 0 || length <= 0 || gravity <= 0 || damping_coefficent < 0)
    {
        return 0;
    }
    if (max_length > length)
    {
        return -1;
    }
    return 1;
}

int 
_min_int(int x, int y)
{
    return x > y ? y : x;
}

int 
_round(double x)
{
    return (int)(x + 0.5);
}

double
_dx(double dx)
{
    return dx;
}

double
_dy(double dy)
{
    return dy;
}

double
_f1(double x, double y, double dx, double tx, double ty, double k, double m, double b, double r)
{
    double L = sqrt(pow(x - tx, 2) + pow(y - ty, 2));
    double s = (x - tx) / L;
    return -(double)k / m * r * s - (float)b / m * dx;
}

double
_f2(double x, double y, double dy, double tx, double ty, double k, double m, double b, double r, double g)
{
    {
        double L = sqrt(pow(x - tx, 2) + pow(y - ty, 2));
        double c = (y - ty) / L;

        return g - (float)k / m * r * c - (float)b / m * dy;
    }
}
// int
// _number_of_files(char *mainDirectoryPath)
// {
//     int subdirectoryCount = 0;
//     DIR *mainDirectory;
//     struct dirent *entry;
//     mainDirectory = opendir(mainDirectoryPath);
//     if (mainDirectory == NULL)
//     {
//         int ret = 0;

//         ret = mkdir(mainDirectoryPath, 0755);

//         if (ret == 0)
//         {
//             printf("Directory created successfully\n");
//             return subdirectoryCount;
//         }
//         else
//         {
//             printf("Unable to create directory %s\n", mainDirectoryPath);
//             return -1;
//         }
//     }
//     while ((entry = readdir(mainDirectory)) != NULL)
//     {
//         if (entry->d_type == DT_DIR)
//         {
//             if (strcmp(entry->d_name, "..") != 0 && strcmp(entry->d_name, ".") != 0)
//             {
//                 subdirectoryCount++;
//             }
//         }
//     }
//     closedir(mainDirectory);
//     return subdirectoryCount;
// }
// int 
// _number_of_files(char *mainDirectoryPath)
// {
//     int subdirectoryCount = 0;
//     DIR *mainDirectory;
//     struct dirent *entry;
//     mainDirectory = opendir(mainDirectoryPath);

//     while ((entry = readdir(mainDirectory)) != NULL)
//     {

//         if (strcmp(entry->d_name, "..") != 0 && strcmp(entry->d_name, ".") != 0)
//         {
//             subdirectoryCount++;
//         }
//     }
//     closedir(mainDirectory);
//     return subdirectoryCount;
// }
// int 
// _directory_create(char *mainDirectoryPath)
// {
//     int subdirectoryCount = 0;
//     DIR *mainDirectory;
//     struct dirent *entry;
//     mainDirectory = opendir(mainDirectoryPath);
//     if (mainDirectory == NULL)
//     {
//         int ret = 0;

//         ret = mkdir(mainDirectoryPath, 0755);

//         if (ret == 0)
//         {
//             printf("Directory created successfully\n");
//             return subdirectoryCount;
//         }
//         else
//         {
//             printf("Unable to create directory %s\n", mainDirectoryPath);
//             return -1;
//         }
//     }
//     closedir(mainDirectory);
//     return 0;
// }
