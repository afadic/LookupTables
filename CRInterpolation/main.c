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
==2284542== I   refs:      4,280,257
==2284542== I1  misses:        1,910
==2284542== LLi misses:        1,853
==2284542== I1  miss rate:      0.04%
==2284542== LLi miss rate:      0.04%
==2284542== 
==2284542== D   refs:      2,006,498  (1,641,871 rd   + 364,627 wr)
==2284542== D1  misses:       29,868  (    7,376 rd   +  22,492 wr)
==2284542== LLd misses:        8,927  (    1,576 rd   +   7,351 wr)
==2284542== D1  miss rate:       1.5% (      0.4%     +     6.2%  )
==2284542== LLd miss rate:       0.4% (      0.1%     +     2.0%  )
==2284542== 
==2284542== LL refs:          31,778  (    9,286 rd   +  22,492 wr)
==2284542== LL misses:        10,780  (    3,429 rd   +   7,351 wr)
==2284542== LL miss rate:        0.2% (      0.1%     +     2.0%  )
==2284542== 
==2284542== Branches:        286,214  (  271,321 cond +  14,893 ind)
==2284542== Mispredicts:       5,640  (    5,321 cond +     319 ind)
==2284542== Mispred rate:        2.0% (      2.0%     +     2.1%   )

**********************************************************************/

#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "time.h" // to compute time
#include "math.h" // math functions
#include "funFile.h" //calls required functions for this file, found in funFile.c
#define MICROSECONDS_IN_SEC 1000000.0

int main(){
    int nDimIn = 8;
    int nDimOut = 1;
    int nVals;
    int i;
    writeTableConfig(nDimIn, nDimOut); //this writes the config file for the table

    double *ptrFirstVal;
	ptrFirstVal = readTableConfig(nDimIn); //pointer to the address of the first entry of table config

    nVals = getNumVals(ptrFirstVal+2, nDimIn); //this gets the number of values of the grid.
    
    printf("Table size %f MB \n",(double) nVals*nDimOut/1000/1024*sizeof(double));

    double *grid;
    grid = writeGrid(ptrFirstVal, nDimIn, nVals); //writes the grid. Not working for more than 6 dimension. Fragmentate memory
    saveGrid(grid,nDimIn,nVals); //save grid to a file. Comment if necessary
    
    //grid = readGrid(nDimIn,nVals); //reads the grid if stored (faster calculations)

    clock_t start = clock(), diff;
    writeTable(nVals,nDimIn,nDimOut,grid); //unnecessary if the table exist
    free(grid);
    diff = clock() - start; int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken write table %d seconds %d milliseconds \n", msec/1000, msec%1000);

    double *ptrTableFirst; ptrTableFirst=readFile(nVals,nDimOut); // only read once. Table stored in memory ready to use

    for(i=0;i<nDimOut;++i){
            printf("table %i first value is: %f \n", i+1, *(ptrTableFirst+nVals*i)); //show table stats
    }

    //testValue = funTransformIn(testValue);
    double sol=0;
    int nIt=1;

    double testValue[nIt][nDimIn];
    int j=0;

    for(i=0;i<nIt;i++) {
        for(j=0;j<nDimIn;j++){
            testValue[i][j] = 1.9;
        }
    }

    clock_t tic = clock();

    for(i=0;i<nIt;i++){
        sol=interpolate(testValue[i],nDimIn,1,ptrFirstVal,ptrTableFirst,nVals);
    }
    clock_t toc = clock();

    free(ptrFirstVal);
    free(ptrTableFirst);

    printf("Average time taken per iteration is %f ms \n", (double)1000*(toc-tic)/CLOCKS_PER_SEC/((double)nIt));
    printf("Solution is %f \n", sol);
    printf("Exact    is %f \n", exactFun(*testValue,nDimIn));

    //Display GNU version
    printf("\n\n\ngcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
    displayAuthor();
    return 0;
}
