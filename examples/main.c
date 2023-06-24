#include "../include/physics.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
int _number_of_files(char *mainDirectoryPath)
{
    int subdirectoryCount = 0;
    DIR *mainDirectory;
    struct dirent *entry;
    mainDirectory = opendir(mainDirectoryPath);

    while ((entry = readdir(mainDirectory)) != NULL)
    {

        if (strcmp(entry->d_name, "..") != 0 && strcmp(entry->d_name, ".") != 0)
        {
            subdirectoryCount++;
        }
    }
    closedir(mainDirectory);
    return subdirectoryCount;
}
int main(void)
{

    printf("%d\n", _number_of_files("."));
}