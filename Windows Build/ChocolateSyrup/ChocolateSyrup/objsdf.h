//
//  objsdf.h
//  ChocolateSyrup
//
//  Created by Karl Li on 4/10/12.
//  Copyright (c) 2012 University of Pennsylvania. All rights reserved.
//

#ifndef ChocolateSyrup_objsdf_h
#define ChocolateSyrup_objsdf_h

#include "makelevelset3.h"
#include "glm/glm.hpp"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

class objsdf
{
private:
    
    glm::vec3 min_box;
    glm::vec3 max_box;
    Array3f phi_grid;
	float scaleFactor;
public:
	std::vector<glm::vec3> vertList;
    std::vector<glm::vec3> faceList;
    objsdf(string filename);
    ~objsdf();
    float getSignedDistanceAtCellObjSpace(int x, int y, int z);
	float getSignedDistanceAtCellWorldSpace(float x, float y, float z);
	glm::vec3 getObjBounds();	//returns object space bounds
	glm::vec3 getWorldBounds();	//returns world space bounds
	glm::vec3 meshCenter;
	float getScaleFactor();	
	float cellsize;
};


#endif
