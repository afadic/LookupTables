/**********************************************************************
Function for building, reading, and evaluating a lookup table
using a Hermite 3rd order spline. Meant to be used as UDF in Fluent.
Author: Anton Fadic
University of Alberta
Chemical Engineering
**********************************************************************/

/*********************************************************************
23/12/2016 First release
06/01/2019 Code maintenance. Checked for memory leaks, readibility and name convention.
Compiles without warnings
gcc -pg main.c funFile.c -o hermite -lm
Compiling with either O1, O2 or O3 flags increases the D1 misses and it does not reduce the execution time

Test results
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes
==99415== 
==99415== HEAP SUMMARY:
==99415==     in use at exit: 0 bytes in 0 blocks
==99415==   total heap usage: 18 allocs, 18 frees, 75,280 bytes allocated
==99415== 
==99415== All heap blocks were freed -- no leaks are possible
==99415== 
==99415== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

valgrind --tool=cachegrind --branch-sim=yes 
Cache misses for 1000 iterations, 4 input dimensions, 1 output dimension
**********************************************************************
==252946== 
==252946== I   refs:      229,410,575
==252946== I1  misses:          2,026
==252946== LLi misses:          1,955
==252946== I1  miss rate:        0.00%
==252946== LLi miss rate:        0.00%
==252946== 
==252946== D   refs:       77,459,911  (64,278,591 rd   + 13,181,320 wr)
==252946== D1  misses:         45,908  (    22,837 rd   +     23,071 wr)
==252946== LLd misses:          9,385  (     1,608 rd   +      7,777 wr)
==252946== D1  miss rate:         0.1% (       0.0%     +        0.2%  )
==252946== LLd miss rate:         0.0% (       0.0%     +        0.1%  )
==252946== 
==252946== LL refs:            47,934  (    24,863 rd   +     23,071 wr)
==252946== LL misses:          11,340  (     3,563 rd   +      7,777 wr)
==252946== LL miss rate:          0.0% (       0.0%     +        0.1%  )
==252946== 
==252946== Branches:       16,293,138  (12,415,142 cond +  3,877,996 ind)
==252946== Mispredicts:       538,478  (   538,195 cond +        283 ind)
==252946== Mispred rate:          3.3% (       4.3%     +        0.0%   )

**********************************************************************/

#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "time.h" // to compute time
#include "math.h" // math functions
#include "funFile.h" //calls required functions for this file, found in funFile.c
#include "krahnertp.h" //calls the kinetic mechanism
#include <omp.h>

#define MICROSECONDS_IN_SEC 1000000.0

#define N_DIM 4
#define N_IT 500

int main(){
    srand(0); //fix seed for reproducibility

    int nDimIn = N_DIM;
    int nDimOut = 1;
    int nVals; 
    int i;

    double n_breaks[N_DIM] = {20,30,30,30};
    double l_bounds[N_DIM] = {800, 1e-4,1e-4,1e-4};
    double u_bounds[N_DIM] = {1400,1e-2,1e-2,1e-2};

    writeTableConfig(nDimIn, n_breaks, l_bounds, u_bounds); //this writes the config file for the table

    double *ptrConfig;
	ptrConfig = readTableConfig(nDimIn); //pointer to the address of the first entry of table config
    nVals = getNumVals(ptrConfig+2, nDimIn);
    printf("Debug: nVals is %d\n", nVals);
   
    printf("Table size %f MB \n",(double) nVals*nDimOut/1000/1024*sizeof(double));

    double *grid;
    grid = writeGrid(ptrConfig, nDimIn, nVals); //writes the grid. 

    //saveGrid(grid,nDimIn,nVals); //save grid to a file. Comment if necessary
    //grid = readGrid(nDimIn,nVals); //reads the grid if stored 

    clock_t start = clock(), diff;
    writeTable(nVals,nDimIn,nDimOut,grid); //unnecessary if the table exist
    free(grid);

    diff = clock() - start; int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("CPU time taken write table %d seconds %d milliseconds \n", msec/1000, msec%1000);

    double *ptrTable; ptrTable=readFile(nVals,nDimOut); // only read once. Table stored in memory ready to use

    i = 0;
    printf("table %i first value is: %f \n", i+1, *(ptrTable+nVals*i)); //show table stats

    int nIt=N_IT;

    double testValue[nIt][nDimIn];
    int j=0;

    double min_val = 0;
    double max_val = 0;

    for(i=0;i<nIt;i++) {
        for(j=0;j<nDimIn;j++){
            if (j==0){
                min_val = 1000;
                max_val = 1200;
            } else {
                min_val = 5e-4;
                max_val = 7e-3;
            }
            testValue[i][j] = min_val + (max_val - min_val)*((double)rand() / RAND_MAX); //generate random test values
        }
    }

    double sol[N_IT];
    double exact[N_IT];

    clock_t tic = clock();
    for(i=0;i<nIt;i++){
        sol[i]=interpolate(testValue[i], ptrConfig, ptrTable);
    }
    clock_t toc = clock();
    printf("Average time taken per inference iteration is %f ms \n", (double)1000*(toc-tic)/CLOCKS_PER_SEC/((double)nIt));

    tic = clock();
    for(i=0;i<nIt;i++){
        exact[i]=exactFun(testValue[i], nDimIn);
        //printf("%f %f ", sol[i], exact[i]);
        //printf("Testing point: %f, %f, %f, %f\n", testValue[i][0], testValue[i][1], testValue[i][2], testValue[i][3]);
    }
    toc = clock();
    printf("Average time taken by sundials is %f ms \n", (double)1000*(toc-tic)/CLOCKS_PER_SEC/((double)nIt));
    
    double error=0;
    error = rmse(sol, exact, nIt);
    printf("RMSE is %f \n", error);
    
    free(ptrConfig);
    free(ptrTable);

    //Display GNU version
    printf("\n\n\ngcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
    return 0;
}
