#include "stdio.h" // in out
#include "math.h" // math functions
#include "stdlib.h" // for dynamic allocation
#include "krahnertp.h" //calls required functions for this file, found in funFile.c


void displayAuthor(){
    printf("Written by Anton Fadic \n");
    printf("University of Alberta \n");
    printf("Version 1.0r \n");
}

double rmse(double *x, double *y, int len){
    //computes the root mean square error between two arrays of length len
    double sum=0.0;
    int i=0;
    for(i=0;i<len;i++){
        sum += (x[i]-y[i])*(x[i]-y[i]);
    }
    return sqrt(sum/len);
}

double exactFun(double *x, int InDim){
    //this is the exact function to compare to
    double outputs=0;
    outputs = r_rates(x);
    //following is debug only
    //printf("NH3: %f\n", outputs);
    //return -x[2]*x[0]*x[2] +x[2]*x[1]+ x[0]*x[1]*x[1]/x[2];
    return outputs;
}

void checkInBoundaries(int *fl, double *t,int nDimIn, double *nBreaks){
    //changes the index *fl and *t if location is out of bounds.
    int i=0;
    for(i=0;i<nDimIn;i++){
        if(*(fl+i)<-1){*(fl+i)=-1;*(t+i)=0;}
        if(*(fl+i)>(*(nBreaks+i))-2){*(fl+i)=*(nBreaks+i)-2;*(t+i)=0;}
    }
}

double CrInterp(double *ptrFll, double t){
//this is the magic 4 point interpolation
double fll, fl, fr, frr;
fll = *ptrFll; fl = *(ptrFll+1); fr = *(ptrFll+2); frr = *(ptrFll+3);
return (1 - 3*t*t + 2*t*t*t )*fl + (3*t*t - 2*t*t*t)*fr + (t - 2*t*t + t*t*t)*(fr-fll)*0.5 + (-t*t + t*t*t)*(frr-fl)*0.5;
        // Lagrange polynomial. Works better for coarser grids, but same performance for finer grids. Looses the C1 continuity property.
        /*
        return (fll * (-1.0/6.0*t*t*t + 1.0/2.0*t*t - 1.0/3.0*t) +
                fl  * ( 1.0/2.0*t*t*t - 1.0*t*t     - 1.0/2.0*t + 1.0) +
                fr  * (-1.0/2.0*t*t*t + 1.0/2.0*t*t + 1.0*t) +
                frr * ( 1.0/6.0*t*t*t - 1.0/6.0*t)
                );
        */
}

double *readFile(int length, int nDimOut){
    FILE * fp;

    double *ptrFile = 0;
	ptrFile = (double*) malloc(nDimOut*length*sizeof(double));
    if(ptrFile==NULL){
            printf("Memory not allocated. Table could not be read. Closing");
            free(ptrFile);
            exit(-1);
    }

    fp = fopen ("Table.bin", "rb");
    if(fp == NULL){
        printf("Error in file opening\n");
        exit(-1);
    }

    printf("Reading table file... \n");
    size_t file_read = fread(ptrFile,sizeof(*ptrFile),nDimOut*length,fp);
    fclose (fp);
    return ptrFile;
}

void writeTable(int length, int nDimIn, int nDimOut, double *grid){
    FILE * fp;
    double *ptrTable = 0;

    ptrTable=(double*) malloc(nDimOut*length*sizeof(double));
    if(ptrTable==NULL){
        printf("memory not allocated. Closing...");
        exit(-1);
    }

    //loop for writing the table. j is number of dimensions and i is entry.
    #pragma omp parallel for schedule(static)
    for(int i=0; i<length; ++i){
        for(int j=0;j<nDimOut;++j){
             *(ptrTable+j*length+i) = exactFun(grid+i*nDimIn, nDimIn);
        }
    }
    fp = fopen ("Table.bin", "wb");
    printf("Writing table file... \n");
    fwrite(ptrTable,sizeof(*ptrTable),nDimOut*length,fp);
    free(ptrTable);
    fclose (fp);
    printf("table written successfully \n");
}

void writeTableConfig(int nDimIn, double *breaks, double *l_bounds, double *u_bounds){
    FILE * fp;
    const int configSize = 3 * nDimIn + 2;
    double config[configSize];
    int i=0;

    //write number of dimensions in
    *config = nDimIn;
    //write the number of dimensions out
    *(config+1) = 1;

    //number of breaks. 0-nDimIn-1
    double n_breaks[nDimIn]; 
    for(i=0; i<nDimIn; i++){
        n_breaks[i] = *(breaks + i);
    }

    for (i = 0; i < nDimIn; i++) {
        *(config + 2 + i) = (double) n_breaks[i];
    }

    double min_vals[nDimIn];
    double max_vals[nDimIn];

    for(i=0; i<nDimIn; i++){
        min_vals[i] = l_bounds[i];
        max_vals[i] = u_bounds[i];
    }

    for(i=0; i<nDimIn; i++){
        *(config+nDimIn+2+i) = min_vals[i]; 
        *(config+2*nDimIn+2+i) = max_vals[i]; //Temperature
    }

    fp = fopen ("tableConfig.bin", "wb");
    if (!fp) {
        perror("Error opening tableConfig.bin for writing");
        exit(-1);
        }
    size_t written = fwrite(config, sizeof(double), configSize, fp);
    if (written != (size_t)configSize) {
        fprintf(stderr, "Error: Only wrote %zu/%d elements\n", written, configSize);
    }
    fclose (fp);
    printf("Config file written successfully \n \n");
}

double *readTableConfig(int nDim){
    FILE * fp;
    int i;

    double *ptr2 = 0;
    ptr2 = (double*) malloc((3*nDim+2)*sizeof(double)); //required number of info: nDimIn, nDimOut, 3*nDim (breaks,min,max)

    if(ptr2==NULL){
            printf("Memory not allocated. Closing readTableConfig");
    exit(0);
    }

    fp = fopen ("tableConfig.bin", "rb");
    printf("Reading config file... \n");
    size_t file_read = fread(ptr2,1,(3*nDim+2)*sizeof(*ptr2),fp);
    fclose (fp);
    printf("Config file read successfully \n \n");
    printf("Number of input dimensions is %i \n", (int) *ptr2);
    printf("Number of output dimensions is %i \n", (int) *(ptr2+1));
    for(i=0;i<nDim;++i){
        printf("Dimension %i, breaks = %i, Min = %f , Max = %f \n", i+1,(int) *(ptr2+i+2), *(ptr2+nDim+i+2), *(ptr2+2*nDim+i+2));
    }
    return ptr2;
}

double *writeGrid(double *ptrFirstVal, int nDimIn, int nVals){
    //This function is in its final version. Please do not touch. 
    int i,j;
    double h[nDimIn]; // breakpoint spacing
    double *nBreaks=0; nBreaks = (ptrFirstVal+2);
    double *min=0; min = (ptrFirstVal + 2 + nDimIn);
    double *max=0; max = (ptrFirstVal + 2 + 2*nDimIn);

    //min = funTransformIn(min);
    //max = funTransformIn(max);

    for(i=0;i<nDimIn;++i){
        *(h+i) = (*(max+i)- *(min+i))/(*(nBreaks+i)-1);
    }
    int sumNBreaks=0; for(i=0;i<nDimIn;++i){
        sumNBreaks += *(nBreaks + i);
    }

    //double rValues[sumNBreaks];
    double *rValues=0; rValues=malloc(sizeof(double)*sumNBreaks);

    /*Objective here is to write the list of possible values as rValues={X1,X2...X_b1,Y1,Y2,...Y_b2}
    with following rule:
    Xi+1=Xi+h; X0= *(min+0);
    Yi+1=Yi+h; Y0= *(min+1); */
    int counter=0; j=0;
    while(j<nDimIn){
        for(i=0;i< *(nBreaks+j);++i){
            *(rValues+counter) = *(min+j) + *(h+j)*i;
            //printf("rValues are %f \n", *(rValues + counter));
            counter++;
        } j++;
    }

    /*The objective here is to write the vectors of *grid values from the rValues vector. This is the loops that fails at nDimIn>6*/
    double *grid=0; 
    grid=(double*) malloc(sizeof(double)*nVals*nDimIn*1); 
    if(grid==NULL) {
        printf("grid memory not allocated. Closing.\n"); exit(1);
    }
    //double grid[nVals*nDimIn];

    int temp=1; j=0; int temp2=0;

    while(j<nDimIn){
        //printf("j is %i \n", j);
        for(i=0;i<nVals;++i){
            *(grid+i*nDimIn+j) = *(rValues + temp2 + (((int) floor(i/temp)) % (int) *(nBreaks+j)) ); //hardest function ever
        }
        temp2 += (int) *(nBreaks+j);
        temp *= (int) *(nBreaks+j);
        j++;
    }

    /* //this code is for visualizing the grid. Uncomment if necessary,
    for(i=0;i<nVals;++i){j=0; while(j<nDimIn){ printf("%f ", *(grid+i*nDimIn+j)); j++;} printf("\n"); } */
    free(rValues);
    return grid;
}

int getNumVals(double *ptrFirstVal, int nDim){
    // *ptrFirstVal is pointing to the number of breaks
    int i;
    int nVals=1;
    double dimsArray[nDim];
    for(i=0;i<nDim;++i){
       dimsArray[i] = *(ptrFirstVal+i);
       nVals = nVals*dimsArray[i];
    }
    //printf("number of grid values per output dimension is: %i \n", nVals);
    return nVals;
}

void saveGrid(double *grid, int nDimIn, int nVals){
    FILE *fp;
    fp = fopen("grid.bin","wb");
    fwrite(grid,sizeof(double),nVals*nDimIn,fp);
    fclose(fp);
}


double interpolate(double *x, double *ptrConfig, double *ptrTableFirst) {
    int nDimIn = (int) *ptrConfig;
    double *nBreaks = ptrConfig + 2;
    double *min = ptrConfig + 2 + nDimIn;
    double *max = ptrConfig + 2 + 2 * nDimIn;
    
    double buffer[256]; //preallocate buffer in stack. May run into stack overflow for larger nDimIn. 256 is set for nDimIn<=4
    //turns out that malloc is much slower than stack allocation, and this function is called many times.
    double h[10];     // Assumes max dimensions <= 10 which is reasonable
    double t[10];
    int fl[10];
    
    // 1. Pre-compute strides for the lookup table to avoid recalculating in loops
    int tableStrides[10];
    tableStrides[0] = 1;
    for(int k=1; k < nDimIn; k++) {
        tableStrides[k] = tableStrides[k-1] * (int)nBreaks[k-1];
    }

    // 2. Calculate coordinates and clamp boundaries
    for (int k = 0; k < nDimIn; k++) {
        h[k] = (max[k] - min[k]) / (nBreaks[k] - 1); 
        
        double u = (x[k] - min[k]) / h[k];
        int idx = (int)floor(u) - 1; 

        // SAFETY: Clamp indices to avoid reading garbage memory (-1)
        if (idx < 0) idx = 0;
        if (idx > (int)nBreaks[k] - 4) idx = (int)nBreaks[k] - 4;

        fl[k] = idx;
        // t is distance from the 2nd point (P1) in the 4-point stencil
        // t should generally be in [0, 1]
        t[k] = (x[k] - min[k]) / h[k] - (double)(fl[k] + 1);
    }

    // 3. Calculate Base Index of the Hypercube
    int baseIndex = 0;
    for(int k=0; k < nDimIn; k++){
        baseIndex += fl[k] * tableStrides[k];
    }

    int refIndex=0; int temp=1; int j=0;
    while(j<nDimIn-1){
        temp *=  ((int) *(nBreaks+j));j++;
        //printf("%i \n", (int) *(fl+j));
        refIndex += temp*(*(fl+j));
    }
    int nVals=0; nVals = getNumVals(ptrConfig + 2, nDimIn);

    refIndex += (*fl); //printf("refIndex is %i \n", refIndex);
    //change the output dimension, include in refIndex
    refIndex = refIndex + (1-1)*nVals; //Read the right table
    int index=0; int coordinate=0; int factor;

    // 4. Extract 4^N points 
    int nPoints = 1 << (2 * nDimIn); // 4^nDimIn
    
    for (int i = 0; i < nPoints; ++i) {
        int tempI = i;
        int currentOffset = 0;
        
        for (int dim = 0; dim < nDimIn; dim++) {
            int local_coord = tempI % 4; // 0, 1, 2, or 3
            tempI /= 4;
            currentOffset += local_coord * tableStrides[dim];
        }
        
        // ptrTableFirst is 1D array. We add the base corner + local offsets
        buffer[i] = ptrTableFirst[baseIndex + currentOffset];
    }

    // 5. Contraction code (Dimension Reduction)
    // Reduce one dimension at a time
    int currentLen = nPoints;
    for (int j = 0; j < nDimIn; j++) {
        currentLen /= 4; // Output size shrinks by 4 each step
        for (int i = 0; i < currentLen; i++) {
            // Interpolate the contiguous block of 4 values using t[j]
            // We overwrite the buffer in-place to save memory
            buffer[i] = CrInterp(buffer + 4 * i, t[j]);
        }
    }

    // The result has contracted to the first element
    return buffer[0];
}
