// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b


// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;



#endif