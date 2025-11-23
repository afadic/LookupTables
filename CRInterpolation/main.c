/**********************************************************************
Function for building, reading, and evaluating a lookup table
using a Hermite 3rd order spline. Meant to be used as UDF in Fluent.
Compiled with GCC on windows.
Author: Anton Fadic
University of Alberta
Chemical Engineering
**********************************************************************/
//#include "udf.h" //UDF macros
#include "stdio.h" // in out
#include "stdlib.h" // for dynamic allocation
#include "time.h" // to compute time
#include "math.h" // math functions
#include "funFile.h" //calls required functions for this file, found in funFile.c
/*********************************************************************
Anton Fadic
23/12/2016 First release
06/01/2019 Code maintenance. Checked for memory leaks, readibility and name convention.
**********************************************************************/

int main(){
    int nDimIn = 2;
    int nDimOut = 1;
    int nVals;
    int i;
    writeTableConfig(nDimIn, nDimOut); //this writes the config file for the table

    double *ptrFirstVal;
	ptrFirstVal = readTableConfig(nDimIn); //pointer to the address of the first entry of table config

    nVals = getNumVals(ptrFirstVal+2, nDimIn); //this gets the number of values of the grid.
    printf("Table size %f MB \n",(double) nVals*nDimOut/1024/1024*sizeof(double));

    double *grid;
    grid = writeGrid(ptrFirstVal, nDimIn, nVals); //writes the grid. Not working for more than 6 dimension. Fragmentate memory
    saveGrid(grid,nDimIn,nVals); //save grid to a file. Comment if necessary

    //grid = readGrid(nDimIn,nVals); //reads the grid if stored (faster calculations)

    clock_t start = clock(), diff;
    writeTable(nVals,nDimIn,nDimOut,grid); //unnecessary if the table exist
    diff = clock() - start; int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken write table %d seconds %d milliseconds \n", msec/1000, msec%1000);

    double *ptrTableFirst; ptrTableFirst=readFile(nVals,nDimOut); // only read once. Table stored in memory ready to use

    for(i=0;i<nDimOut;++i){
            printf("table %i first value is: %f \n", i+1, *(ptrTableFirst+nVals*i)); //show table stats
    }

    double testValue[2];

    testValue[0] = 1.9;  //temp
    testValue[1] = 1.9;  //yNH3
    // *(testValue + 2) = 0.05;  //yNO
    // *(testValue + 3) = 0.2;  //yO2

    //testValue = funTransformIn(testValue);
    double *sol=0;
    int nIt=1;
    clock_t tic = clock();
    for(i=0;i<nIt;i++){
        //interpolate(testValue,nDimIn,1,ptrFirstVal,ptrTableFirst,nVals); // receives the pointer of the stored interpolated data. Order is: nDimIn,nDimOut,numTestVal
        sol=interpolate(testValue,nDimIn,1,ptrFirstVal,ptrTableFirst,nVals);
    }
    clock_t toc = clock();
    printf("Average time taken per iteration is %f ms \n", (double)1000*(toc-tic)/CLOCKS_PER_SEC/((double)nIt));
    printf("Solution is %f \n", *sol);
    printf("Exact    is %f \n", exactFun(testValue,1,1));

    //Display GNU version
    printf("\n\n\ngcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
    displayAuthor();
    return 0;
}
