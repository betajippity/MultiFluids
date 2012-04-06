// Based on sample code from Paul Bourke

#ifndef MARCHING_CUBES_H_
#define MARCHING_CUBES_H_

#include <vector>
#include "array3.h"
#include "glm/glm.hpp"
#include <sstream>

float fGetOffset(float fValue1, float fValue2, float fValueDesired);

void MarchingCubes(const Array3f& liquid, float dx, int frameNum, bool outputOBJ);
void MarchCube(const Array3f& liquid, float fX, float fY, float fZ, float dx);

#endif  // MARCHING_CUBES_H_