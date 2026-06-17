#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"
#include <stdlib.h>

int main(int argc, char* argv[]){

    //For timing the code
    clock_t start_time = clock();
    clock_t end_time;
    double cpu_time_used;

    time_t start_wall = time(NULL);

  

    */PUT REST OF MAIN CODE HERE/*


      

    //End bit of timing code
    //Put as final bit of main function
    end_time = clock();
    time_t end_wall = time(NULL);

    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nTotal (CPU) simulation time: %f seconds\n", cpu_time_used);

    printf("Wall clock time: %ld seconds\n", end_wall - start_wall);

  return 0;
}
