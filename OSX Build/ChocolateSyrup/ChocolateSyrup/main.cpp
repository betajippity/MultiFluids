#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>
#include <stdio.h>
#include <cmath>
#include <string>

#include "camera.h"
#include "fps.h"
#include "constants.h"

#include "fluidsim.h"
#include "array3_utils.h"
#include "glm/glm.hpp"
#include "open_gl_headers.h"
#include "basic_math.h"

using namespace std;

//-------------------------------------------------------------
//----------------------SET UP STUFF---------------------------
//-------------------------------------------------------------

string outpath;

int grid_resolution = 10;
float timestep = 0.01f;
int frame = 0;

float grid_width = 1;

bool textOutput = false;

FluidSim sim;

Camera theCamera;
mmc::FpsTracker theFpsTracker;

// UI Helpers
int lastX = 0, lastY = 0;
int theMenu = 0;
int theButtonState = 0;
int theModifierState = 0;
bool isRunning = true;

int savedWidth = 0;
int savedHeight = 0;

float ground_height = 0.2f;

//-------------------------------------------------------------
//-------------------BOUNDARY CONDITIONS-----------------------
//-------------------------------------------------------------

float ground_phi(const glm::vec3& position, float height) {
  return height - position.y;
}

float rad0 = 0.35f;

float box_phi(const glm::vec3& position, glm::vec3& centre, glm::vec3 halfDim) {
  float distX = abs(position.x - centre.x) - halfDim.x;
  float distY = abs(position.y - centre.y) - halfDim.y;
  float distZ = abs(position.z - centre.z) - halfDim.z;

    return max(distX, max(distY, distZ));
}

float sphere_phi(const glm::vec3& position, const glm::vec3& centre, float radius) {
    return (glm::length(position-centre) - radius);
}

float sphereInBox_phi(const glm::vec3& position, glm::vec3& centre) {
  float box = box_phi(position, centre, glm::vec3(rad0, rad0, rad0));

  float sphere = -sphere_phi(position, centre, 0.30f);

  return max(box, sphere);
}

float obj_phi(const glm::vec3& position){
    std::cout << position[0] << " " << position[1] << " " << position[2] << std::endl;
}


float boundary_phi(const glm::vec3& position) {
    
    glm::vec3 c0(0.5f,0.5f,0.5f);
    
    //Using OBJ SDF constraint
    //return obj_phi(position);
    
    //Using the sphere in box constraint
    return -sphereInBox_phi(position, c0);

    // Using the ground constraint.
    //return -ground_phi(position, ground_height);

    // Using the box constraint.
    //return -box_phi(position, c0, glm::vec3(rad0, rad0, rad0));

    // Using the sphere constraint. 
    //return sphere_phi(position, c0, rad0);
}

float liquid_phi(const glm::vec3& position);
float liquid_phi(const glm::vec3& position) {
  
   return sphere_phi(position/grid_width, glm::vec3(0.55f, 0.8f, 0.4f), 0.23f);
   
}

//-------------------------------------------------------------
//----------------------OUTPUT STUFF---------------------------
//-------------------------------------------------------------

void export_particles(string path, int frame, const std::vector<particle*>& particles, float radius) {
    //Write the output
    
    std::stringstream strout;
    strout << path << "particles_" << frame << ".txt";
    string filepath = strout.str();
    
    ofstream outfile(filepath.c_str());
    //write vertex count and particle radius
    outfile << particles.size() << " " << radius << std::endl;
    //write vertices
    for(unsigned int i = 0; i < particles.size(); ++i)
        outfile << particles[i]->position[0] << " " << particles[i]->position[1] << " " << particles[i]->position[2] << std::endl;
    outfile.close();
}

//-------------------------------------------------------------
//----------------------OPENGL STUFF---------------------------
//-------------------------------------------------------------

void initCamera()
{
   double w = theDim[0]*theCellSize;   
   double h = theDim[1]*theCellSize;   
   double d = theDim[2]*theCellSize;   
   double angle = 0.5*theCamera.dfltVfov*BasicMath::PI/180.0;
   double dist;
   if (w > h) dist = w*0.5/std::tan(angle);  // aspect is 1, so i can do this
   else dist = h*0.5/std::tan(angle);
  // theCamera.dfltEye = glm::vec3(w*0.5, h*0.5, -(dist+d*0.5));

   theCamera.dfltEye = glm::vec3(0,0,1);
   //theCamera.dfltLook = glm::vec3(w*0.5, h*0.5, 0.0);
   theCamera.dfltLook = glm::vec3(0, 0, 0.0);
   theCamera.reset();
}

void onMouseMotionCb(int x, int y)
{
   int deltaX = lastX - x;
   int deltaY = lastY - y;
   bool moveLeftRight = abs(deltaX) > abs(deltaY);
   bool moveUpDown = !moveLeftRight;

   if (theButtonState == GLUT_LEFT_BUTTON)  // Rotate
   {
      if (moveLeftRight && deltaX > 0) theCamera.orbitLeft(deltaX);
      else if (moveLeftRight && deltaX < 0) theCamera.orbitRight(-deltaX);
      else if (moveUpDown && deltaY > 0) theCamera.orbitUp(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.orbitDown(-deltaY);
   }
   else if (theButtonState == GLUT_MIDDLE_BUTTON) // Zoom
   {
   //   if (moveUpDown && deltaY > 0) theCamera.moveForward(deltaY);
    //  else if (moveUpDown && deltaY < 0) theCamera.moveBack(-deltaY);
   }    

   if (theModifierState & GLUT_ACTIVE_ALT) // camera move
   {
      if (theButtonState == GLUT_RIGHT_BUTTON) // Pan
      {
         if (moveLeftRight && deltaX > 0) theCamera.moveLeft(deltaX);
         else if (moveLeftRight && deltaX < 0) theCamera.moveRight(-deltaX);
         else if (moveUpDown && deltaY > 0) theCamera.moveUp(deltaY);
         else if (moveUpDown && deltaY < 0) theCamera.moveDown(-deltaY);
      }   
   }
 
   lastX = x;
   lastY = y;
   glutPostRedisplay();
}

void onMouseCb(int button, int state, int x, int y)
{
   theButtonState = button;
   theModifierState = glutGetModifiers();
   lastX = x;
   lastY = y;

   //glutSetMenu(theMenu);
   //if (theModifierState & GLUT_ACTIVE_ALT)
   //{
    //  glutDetachMenu(GLUT_RIGHT_BUTTON);
   //}
   //else
  // {
   //   glutAttachMenu(GLUT_RIGHT_BUTTON);
   //}

   onMouseMotionCb(x, y);
}

void onKeyboardCb(unsigned char key, int x, int y)
{
   if (key == ' ') theCamera.reset();
   //else if (key == '0') MACGrid::theRenderMode = MACGrid::CUBES;
   //else if (key == '1') MACGrid::theRenderMode = MACGrid::SHEETS;
   //else if (key == 'v') MACGrid::theDisplayVel = !MACGrid::theDisplayVel;
   else if (key == 'r') sim.setRecording(!sim.isRecording(), savedWidth, savedHeight);
   else if (key == '>') isRunning = true;
   else if (key == '=') isRunning = false;
   else if (key == '<') sim.reset(grid_width, grid_resolution, grid_resolution, grid_resolution, liquid_phi);
   else if (key == 'q') sim.set_liquid(liquid_phi, glm::vec3(1,0,0));
   else if (key == 'p') sim.set_liquid(liquid_phi, glm::vec3(0,1,0));
   else if (key == 'l') sim.set_liquid(liquid_phi, glm::vec3(0,0,1));
   else if (key == 'w') sim.setTransparentRender(!sim.isTransparentRender());
   else if (key == 'e') sim.setVerbose(!sim.isVerbose());
   else if (key == 't') textOutput = !textOutput;
   else if (key == 'v') sim.outputOBJ = !sim.outputOBJ;
   else if (key == 27) exit(0); // ESC Key
   glutPostRedisplay();
}

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   //case -6: theSmokeSim.reset(); break;
   default: onKeyboardCb(value, 0, 0); break;
   }
}

void onKeyboardSpecialCb(int key, int x, int y)
{
}

void onTimerCb(int value)
{
  if(sim.isVerbose()){
    printf("--------------------\nFrame %d\n", sim.getTotalFrames());
    printf("Simulating liquid\n");
   }
   if (isRunning){
     sim.advance(timestep);
   }
   if(textOutput){
     export_particles(outpath, sim.getTotalFrames(), sim.particles, sim.particle_radius);
     if(sim.isVerbose()){
      printf("Exporting particle data\n");
     }
   }

   glutTimerFunc(theMillisecondsPerFrame, onTimerCb, 0);
   glutPostRedisplay();
}

void onResizeCb(int width, int height)
{
  // Save the width and height:
  savedWidth = width;
  savedHeight = height;
  
   // Update viewport
   glViewport(0, 0, width, height);

   // Update camera projection's aspect ratio
   float vfov, aspect, zNear, zFar;
   theCamera.getProjection(&vfov, &aspect, &zNear, &zFar);
   theCamera.setProjection(vfov, ((GLfloat) width)/height, zNear, zFar);
}

void drawOverlay()
{
  // Draw Overlay
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushAttrib(GL_LIGHTING_BIT);
     glDisable(GL_LIGHTING);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);

     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
     glRasterPos2f(0.01, 0.01); 
     
     char info[1024];
     sprintf(info, "CHOCOLATE SYRUP | Framerate: %3.1f  |  Frame: %u  |  %s", 
         theFpsTracker.fpsAverage(), sim.getTotalFrames(),//,
         sim.isRecording()? "Recording..." : "");
 
     for (unsigned int i = 0; i < strlen(info); i++)
     {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
     }
  glPopAttrib();
}

void drawAxes()
{
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);

    glLineWidth(2.0); 
    glBegin(GL_LINES);
      glColor3f(1.0, 0.0, 0.0);
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(1.0, 0.0, 0.0);

      glColor3f(0.0, 1.0, 0.0);
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(0.0, 1.0, 0.0);

      glColor3f(0.0, 0.0, 1.0);
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(0.0, 0.0, 1.0);
    glEnd();
  glPopAttrib();
}

void onDrawCb()
{
  // Keep track of time
  theFpsTracker.timestamp();

  // Draw Scene and overlay
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  theCamera.draw();
  drawAxes();
  sim.draw();
  drawOverlay();
  glutSwapBuffers();
}

void init(void)
{
    initCamera();
    glClearColor(.4, .4, .4, 1.0);

    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);
    glDisable(GL_LIGHTING);
    glCullFace(GL_BACK);
}

int main(int argc, char **argv)
{

  if(argc!=2){
      cerr << "The first parameter should be the folder to write the output liquid meshes into. (eg. c:\\output\\)" << endl;
      return 1;
   }

   string output(argv[1]);

   outpath = output;

   sim.initialize(grid_width, grid_resolution, grid_resolution, grid_resolution);
   sim.set_boundary(boundary_phi);
   sim.set_liquid(liquid_phi, glm::vec3(0,0,1));

   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
   glutInitWindowSize(1280, 720);
   glutInitWindowPosition(100, 100);
   glutCreateWindow("Chocolate Syrup - CIS563");
   glutDisplayFunc(onDrawCb);
   glutKeyboardFunc(onKeyboardCb);
   glutSpecialFunc(onKeyboardSpecialCb);
   glutMouseFunc(onMouseCb);
   glutMotionFunc(onMouseMotionCb); 
   glutTimerFunc(theMillisecondsPerFrame, onTimerCb, 0);
   glutReshapeFunc(onResizeCb);

    /*int viewMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Toggle velocities\t'v'", 'v');
    glutAddMenuEntry("Render density as cubes\t'0'", '0');
    glutAddMenuEntry("Render density as sheets\t'1'", '1');
  glutAddMenuEntry("Render color from density\t'd'", '2');
  glutAddMenuEntry("Render color from temperature\t't'", '3');

    theMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Start\t'>'", '>');
    glutAddMenuEntry("Pause\t'='", '=');
    glutAddMenuEntry("Reset\t'<'", '<');
    glutAddMenuEntry("Reset camera\t' '", ' ');
    glutAddMenuEntry("Record\t'r'", 'r');
    glutAddSubMenu("Display", viewMenu);
    glutAddMenuEntry("_________________", -1);
    glutAddMenuEntry("Exit", 27);
    glutAttachMenu(GLUT_RIGHT_BUTTON);*/

    init();

    glutMainLoop();
    return 0;   
}

