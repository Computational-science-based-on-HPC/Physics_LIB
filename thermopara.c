#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "thermopara.h"
#include "thermoutils.h"
#define ll long long
#define M_PI 3.14159265358979323846264338327

// #ifdef __linux__
#include "mpi.h"
#include "omp.h"
// #endif

double _get_value(struct TimeParam* time, struct SpaceParam* space, double x, double t, int precision)
{
    double sum = 0.0, exponential, spaceXTerm, coeff;
    double x_real = x * (space->delta_x);
    double t_real = t * time->delta_t;

    for (int k = 0; k < precision; k++)
    {
        exponential = exp(-3 * pow(2 * (k + 1), 2) * (M_PI * M_PI * t_real) / 4);
        spaceXTerm = sin((double)(2 * k + 1) * M_PI * x_real / 2);
        coeff = 1 / (2 * k + 1);
        sum += coeff * exponential * spaceXTerm;
    }

    sum *= 200 / M_PI;
    return sum;
}

int _simulate_heat_transfer_1D_MPI(struct TimeParam* time_param, struct SpaceParam* space_param, int precision){
    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    fptr1 = fopen("1.txt", "w");
    fptr2 = fopen("2.txt", "w");
    fptr3 = fopen("3.txt", "w");
    int my_rank;     // rank of process
    int processesNo; // number of process

    int checkRem;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesNo);

    ll numTimePoint;
    _cal_num_time(time_param, &numTimePoint);

    ll numSpacePoint;
    _cal_num_space(space_param, &numSpacePoint);

    ll numTimePointPerProcess = numTimePoint / (processesNo - 1); // number of time points per process
    ll numTimePointRemProcess = numTimePoint % (processesNo - 1); // number of time points for last process

    MPI_Bcast(&numTimePointPerProcess, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSpacePoint, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&space_param, sizeof(space_param), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_param, sizeof(time_param), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        ll startIndex = 0;
        int i;
        for (i = 1; i < processesNo; i++)
        {
            int flag = 0;
            if (numTimePointRemProcess > 0)
            {
                flag = 1;
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess + 1;
                numTimePointRemProcess -= 1;
            }
            else
            {
                MPI_Send(&flag, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&startIndex, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD);
                startIndex += numTimePointPerProcess;
            }
        }
    }
    else
    {
        ll startIndex;
        ll i;

        MPI_Recv(&checkRem, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        ll size;
        if (checkRem > 0)
        {
            size = numTimePointPerProcess + 1;
        }
        else
        {
            size = numTimePointPerProcess;
        }

        MPI_Recv(&startIndex, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);

        printf("The value of startIndex is: %lld\n", startIndex);

        i = startIndex;
        ll endIndex = startIndex + size;

        for (; i < endIndex; ++i)
        {
            for (ll x = 0; x < numSpacePoint; ++x)
            {
                if (my_rank == 1)
                {
                    fprintf(fptr1, "%f ", getValue(&time_param, &space_param, x, i));
                }
                else if (my_rank == 2)
                {
                    fprintf(fptr2, "%f ", getValue(&time_param, &space_param, x, i));
                }
                else if (my_rank == 3)
                {
                    fprintf(fptr3, "%f ", getValue(&time_param, &space_param, x, i));
                }
            }
            if (my_rank == 1)
            {
                fprintf(fptr1, "\n");
            }
            else if (my_rank == 2)
            {
                fprintf(fptr2, "\n");
            }
            else if (my_rank == 3)
            {
                fprintf(fptr3, "\n");
            }
        }
    }
    fclose(fptr1);
    fclose(fptr2);
    fclose(fptr3);
    MPI_Finalize();
    return 0;
}
