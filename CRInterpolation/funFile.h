
double exactFun(double *inVal, int outDim, int InDim); //exact function to be approximated

double CrInterp(double *ptrFll, double t); //interpolates

double *readFile(int length, int nDimOut); //read the table

void writeTable(int length, int nDimOut, int nDimIn, double *grid); // writes the table

void writeTableConfig(int nDim, int nDimOut); //write config file of table.

double *readTableConfig(int nDim); //read the table config

double *writeGrid(double *ptrFirstVal, int nDimIn, int nVals); // writes the grid

int getNumVals(double *ptrFirstVal, int nDim); // get the number of values from the number of dimentions

void saveGrid(double *grid, int nDimIn, int nVals); // stores the grid in Disk

void *interpolate(double *x, int nDimIn, int nDimOut, double *ptrFirstVal, double *ptrTableFirst, int numTestVal); //the most important function in the world

double *funTransformIn(double *inVal);

double *invFunTransformIn(double *inVal);

void checkInBoundaries(int *fl, double *t,int nDimIn);

int *getPosition(int *coordinate, int *nBreaks, int nDimIn);

void displayAuthor();
