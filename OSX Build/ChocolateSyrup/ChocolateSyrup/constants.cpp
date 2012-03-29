// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {11, 11, 1};
#else
//const int theDim[3] = {12, 12, 4};
const int theDim[3] = {51, 51, 1};
#endif

const double theCellSize = 1.0; //0.5;

