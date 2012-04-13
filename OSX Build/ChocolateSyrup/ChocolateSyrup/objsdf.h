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

using namespace std;

class objsdf
{
private:
    std::vector<glm::vec3> vertList;
    std::vector<glm::vec3> faceList;
    glm::vec3 min_box;
    glm::vec3 max_box;
    Array3f phi_grid;
public:
    objsdf(string filename);
    ~objsdf();
    float getSignedDistanceAtCell(int x, int y, int z);
};


#endif
