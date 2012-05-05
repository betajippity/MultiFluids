#include "fluidsim.h"
#include "levelset_util.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "glm/glm.hpp"
#include <omp.h>
#include "stb_image_write.h"
#include "marching_cubes.h"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

struct particle;

void extrapolate(Array3f& grid, Array3c& valid);

FluidSim::FluidSim(){
	//mesh = new objsdf("C:/Users/Karl/Dropbox/Documents/Projects/SDFGen/cow.obj");
	//mesh = new objsdf("C:/Users/Karl/Dropbox/Documents/Projects/SDFGen/bunny_watertight.obj");
}
FluidSim::~FluidSim(){
    
}


int FluidSim::getTotalFrames(){
	return mTotalFrameNum;
}

bool FluidSim::isRecording()
{
   return mRecordEnabled;
}

void FluidSim::setRecording(bool on, int width, int height)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
	
	recordWidth = width;
	recordHeight = height;
}

void FluidSim::grabScreen()
{
	
	if (mFrameNum > 9999) exit(0);
	
	// Save an image:

	unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

	for (int i=0; i<recordHeight; i++) 
	{
		glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
			bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
	}
	 
	char anim_filename[2048];
	//sprintf_s(anim_filename, 2048, "../output/images/fluid_%04d.png", mFrameNum); 
	sprintf(anim_filename, "../output/images/fluid_%04d.png", mFrameNum); 
    
	stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
	
	delete [] bitmapData;
	
	mFrameNum++;
	 
}

void FluidSim::initialize(float width, int ni_, int nj_, int nk_) {
	frameNum = 0;
	outputOBJ = false;
   ni = ni_;
   nj = nj_;
   nk = nk_;
   dx = width / (float)ni;
    grid_width = width; 
    
   u.resize(ni+1,nj,nk); temp_u.resize(ni+1,nj,nk); u_weights.resize(ni+1,nj,nk); u_valid.resize(ni+1,nj,nk);
   u_vol.resize(ni+1, nj, nk);
   v.resize(ni,nj+1,nk); temp_v.resize(ni,nj+1,nk); v_weights.resize(ni,nj+1,nk); v_valid.resize(ni,nj+1,nk);
   v_vol.resize(ni, nj+1, nk);
   w.resize(ni,nj,nk+1); temp_w.resize(ni,nj,nk+1); w_weights.resize(ni,nj,nk+1); w_valid.resize(ni,nj,nk+1);
   w_vol.resize(ni, nj, nk+1);

   c_vol.resize(ni, nj, nk);
   n_vol.resize(ni+1, nj+1, nk+1);

   viscosity.resize(ni, nj, nk);
   viscosity.assign(0.0f);

   particle_radius = (float)(dx * 1.01*sqrt(3.0)/2.0); 
   //make the particles large enough so they always appear on the grid

   u.set_zero();
   v.set_zero();
   w.set_zero();
   nodal_solid_phi.resize(ni+1,nj+1,nk+1);
   valid.resize(ni+1, nj+1, nk+1);
   old_valid.resize(ni+1, nj+1, nk+1);
   liquid_phi.resize(ni,nj,nk);
   mRecordEnabled=false;
   transparentRender=true;
}

void FluidSim::setVerbose(bool set){
	verbose = set;
}

bool FluidSim::isVerbose(){
	return verbose;
}

void FluidSim::setTransparentRender(bool set){
	transparentRender = set;
}

bool FluidSim::isTransparentRender(){
	return transparentRender;
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim::set_boundary(float (*phi)(const glm::vec3&)) {

   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      glm::vec3 pos(i*dx,j*dx,k*dx);
      nodal_solid_phi(i,j,k) = phi(pos);
   }

}

void FluidSim::set_liquid(float (*phi)(const glm::vec3&), glm::vec3 color, double viscosityCoef) {
   //surface.reset_phi(phi, dx, Vec3f(0.5f*dx,0.5f*dx,0.5f*dx), ni, nj, nk);
   
   //initialize particles
   int seed = 0;
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      glm::vec3 pos(i*dx,j*dx,k*dx);
      float a = randhashf(seed++); float b = randhashf(seed++); float c = randhashf(seed++);
      pos += dx * glm::vec3(a,b,c);

      if(phi(pos) <= -particle_radius) {
         glm::vec3 posdividedbydx = pos/dx;
         
         float solid_phi =  interpolate_value<glm::vec3>(posdividedbydx, nodal_solid_phi);
         if(solid_phi >= 0){
			particle* pt = new particle();
			pt->position = pos;
			pt->color = color;
			pt->viscosity = viscosityCoef;
			particles.push_back(pt);		 
			//particles.push_back(pos);
			//colors.push_back(color);
		 }
      }
   }
}

void FluidSim::reset(float width, int ni_, int nj_, int nk_, float (*phi)(const glm::vec3&))
{
    initialize(width, ni_, nj_, nk_);
    particles.clear();
    set_liquid(phi, glm::vec3(0,0,1), 0.0f);
	mTotalFrameNum = 0;
}

//The main fluid simulation step
void FluidSim::advance(float dt) {
   float t = 0;

   while(t < dt) {
      float substep = cfl();   
      if(t + substep > dt)
         substep = dt - t;
      //printf("Taking substep of size %f (to %0.3f%% of the frame)\n", substep, 100 * (t+substep)/dt);
      
      //printf(" Surface (particle) advection\n");
      advect_particles(substep);

      //printf(" Velocity advection\n");
      //Advance the velocity
      advect(substep);
      add_force(substep);

	  //apply_viscosity(substep);

      //printf(" Pressure projection\n");
      project(substep); 
       
      //Pressure projection only produces valid velocities in faces with non-zero associated face area.
      //Because the advection step may interpolate from these invalid faces, 
      //we must extrapolate velocities from the fluid domain into these invalid faces.
      //printf(" Extrapolation\n");
      extrapolate(u, u_valid);
      extrapolate(v, v_valid);
      extrapolate(w, w_valid);
    
      //For extrapolated velocities, replace the normal component with
      //that of the object.
      //printf(" Constrain boundary velocities\n");
      constrain_velocity();

      t+=substep;
   }

   mTotalFrameNum++;
}


float FluidSim::cfl() {

   float maxvel = 0;
   for(unsigned int i = 0; i < u.a.size(); ++i)
      maxvel = max(maxvel, (float)fabs(u.a[i]));
   for(unsigned int i = 0; i < v.a.size(); ++i)
      maxvel = max(maxvel, (float)fabs(v.a[i]));
   for(unsigned int i = 0; i < w.a.size(); ++i)
      maxvel = max(maxvel, (float)fabs(w.a[i]));
   
   return dx / maxvel;
}

void FluidSim::add_particle(const glm::vec3& pos, const glm::vec3& color) {
	particle* pt = new particle();
	pt->position = pos;
	pt->color = color;
	particles.push_back(pt);
	//particles.push_back(pos);
	//colors.push_back(color);
}

void FluidSim::add_force(float dt) {

   //gravity
   for(int k = 0;k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      v(i,j,k) -= 9.81f * dt;
   }

}



//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim::constrain_velocity() {
   temp_u = u;
   temp_v = v;
   temp_w = w;

   //(At lower grid resolutions, the normal estimate from the signed
   //distance function can be poor, so it doesn't work quite as well.
   //An exact normal would do better if we had it for the geometry.)

   //constrain u
	#pragma omp parallel for
   for(int k = 0; k < u.nk;++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      if(u_weights(i,j,k) == 0) {
         //apply constraint
          glm::vec3 pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
          glm::vec3 vel = get_velocity(pos);
          glm::vec3 normal(0,0,0);
          interpolate_gradient<glm::vec3>(normal, pos/dx, nodal_solid_phi); 
          normal = glm::normalize(normal);
          float perp_component = glm::dot(vel, normal);
          vel -= perp_component*normal;
          temp_u(i,j,k) = vel[0];
      }
   }

   //constrain v
   #pragma omp parallel for
   for(int k = 0; k < v.nk;++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      if(v_weights(i,j,k) == 0) {
          //apply constraint
          glm::vec3 pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
          glm::vec3 vel = get_velocity(pos);
          glm::vec3 normal(0,0,0);
          interpolate_gradient<glm::vec3>(normal, pos/dx, nodal_solid_phi); 
          normal = glm::normalize(normal);
          float perp_component = glm::dot(vel, normal);
          vel -= perp_component*normal;
          temp_v(i,j,k) = vel[1];
      }
   }

   //constrain w
   #pragma omp parallel for
   for(int k = 0; k < w.nk;++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      if(w_weights(i,j,k) == 0) {
         //apply constraint
         glm::vec3 pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
         glm::vec3 vel = get_velocity(pos);
         glm::vec3 normal(0,0,0);
         interpolate_gradient<glm::vec3>(normal, pos/dx, nodal_solid_phi); 
         normal = glm::normalize(normal);
         float perp_component = glm::dot(vel, normal);
         vel -= perp_component*normal;
         temp_w(i,j,k) = vel[2];
      }
   }

   //update
   u = temp_u;
   v = temp_v;
   w = temp_w;

}

void FluidSim::advect_particles(float dt) { 
   for(unsigned int p = 0; p < particles.size(); ++p) {
      particles[p]->position = trace_rk2(particles[p]->position, dt);
   
      //check boundaries and project exterior particles back in
      float phi_val = interpolate_value<glm::vec3>(particles[p]->position/dx, nodal_solid_phi); 
      if(phi_val < 0) {
         glm::vec3 grad;
         interpolate_gradient<glm::vec3>(grad, particles[p]->position/dx, nodal_solid_phi);
          if(glm::length(grad) > 0)
            grad = glm::normalize(grad);
         particles[p]->position -= phi_val * grad;
      }
   }
   

}

//Basic first order semi-Lagrangian advection of velocities
void FluidSim::advect(float dt) {

   temp_u.assign(0);
   temp_v.assign(0);
   temp_w.assign(0);

   #pragma omp parallel for
   //semi-Lagrangian advection on u-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
       glm::vec3 pos(i*dx, (j+0.5f)*dx, (k+0.5f)*dx);
       pos = trace_rk2(pos, -dt);
       temp_u(i,j,k) = get_velocity(pos)[0];  
   }

   #pragma omp parallel for
   //semi-Lagrangian advection on v-component of velocity
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
       glm::vec3 pos((i+0.5f)*dx, j*dx, (k+0.5f)*dx);
       pos = trace_rk2(pos, -dt);
       temp_v(i,j,k) = get_velocity(pos)[1];
   }

   #pragma omp parallel for
   //semi-Lagrangian advection on w-component of velocity
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
       glm::vec3 pos((i+0.5f)*dx, (j+0.5f)*dx, k*dx);
       pos = trace_rk2(pos, -dt);
       temp_w(i,j,k) = get_velocity(pos)[2];
   }

   //move update velocities into u/v vectors
   u = temp_u;
   v = temp_v;
   w = temp_w;
}

void FluidSim::compute_phi() {
   
   //grab from particles
   liquid_phi.assign(3*dx);
   viscosity.assign(0.0f);
   Array3f viscousCounts;
   viscousCounts.resize(ni,nj,nk);
   viscousCounts.assign(0);
   for(unsigned int p = 0; p < particles.size(); ++p) {
      glm::vec3 cell_ind(particles[p]->position / dx);
      for(int k = max((float)0,(float)cell_ind[2] - 1); k <= min((float)cell_ind[2]+1,(float)nk-1); ++k) {
         for(int j = max((float)0,(float)cell_ind[1] - 1); j <= min((float)cell_ind[1]+1,(float)nj-1); ++j) {
            for(int i = max((float)0,(float)cell_ind[0] - 1); i <= min((float)cell_ind[0]+1,(float)ni-1); ++i) {
				viscosity(i,j,k) = (float)viscosity(i,j,k) + (float)particles[p]->viscosity;
				viscousCounts(i,j,k) = viscousCounts(i,j,k) + 1;
               glm::vec3 sample_pos((i+0.5f)*dx, (j+0.5f)*dx,(k+0.5f)*dx);
                float test_val = glm::length(sample_pos-particles[p]->position) - particle_radius;
               if(test_val < liquid_phi(i,j,k))
                  liquid_phi(i,j,k) = test_val;
            }
         }
      }
   }

    for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) {
         for(int i = 0; i < ni; ++i) {
			 if(viscousCounts(i,j,k) != 0) viscosity(i,j,k) = viscosity(i,j,k) / viscousCounts(i,j,k);
			 else viscosity(i,j,k) = 0;
		 }
	  }
   }
   
   //extend phi slightly into solids (this is a simple, naive approach, but works reasonably well)
   Array3f phi_temp = liquid_phi;
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) {
      for(int j = 0; j < nj; ++j) {
         for(int i = 0; i < ni; ++i) {
            if(liquid_phi(i,j,k) < 0.5*dx) {
               float solid_phi_val = 0.125f*(nodal_solid_phi(i,j,k) + nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i+1,j+1,k)
                  + nodal_solid_phi(i,j,k+1) + nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i+1,j+1,k+1));
               if(solid_phi_val < 0)
                  phi_temp(i,j,k) = -0.5f*dx;
            }
         }
      }
   }
   liquid_phi = phi_temp;
   

}

void FluidSim::apply_viscosity(float dt) {
	compute_phi();

	compute_viscosity_weights(dt);

	//std::cout << "About to solve viscosity" << std::endl;
	solve_viscosity(dt);
	//std::cout << "Done solving viscosity" <<std::endl;
}

int FluidSim::u_ind(int i, int j, int k) {
	return i + j*(ni+1) + (ni+1)*nj*k;
}

int FluidSim::v_ind(int i, int j, int k) {
	return i + j*ni + ni*(nj+1)*k + (ni+1)*nj*nk;
}

int FluidSim::w_ind(int i, int j, int k) {
	return i + j*ni + k*ni*nj + (ni+1)*nj*nk + ni*(nj+1)*nk;
}

void FluidSim::solve_viscosity(float dt) {
   int ni = liquid_phi.ni;
   int nj = liquid_phi.nj;
   int nk = liquid_phi.nk;

   float u_obj = 0;
   float v_obj = 0;
   float w_obj = 0;

   Array3c u_state(ni+1,nj,nk,(const char&)0);
   Array3c v_state(ni,nj+1,nk,(const char&)0);
   Array3c w_state(ni,nj,nk+1,(const char&)0);
   const int SOLID = 1;
   const int FLUID = 0;

   //printf("Determining states\n");

   // Determine the solid boundaries for the x-faces.
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) {
	   for(int j = 0; j < nj; ++j) {
		  for(int i = 0; i < ni+1; ++i) {
			 if(i < 1 || i >= ni || (nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i,j,k)
										+nodal_solid_phi(i,j+1,k+1) + nodal_solid_phi(i,j,k+1))/4.0 <= 0)
				u_state(i,j,k) = SOLID;
			 else
				u_state(i,j,k) = FLUID;
		  }
	   }
   }

   // Determine the solid boundaries for the y-faces. 
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) {
	   for(int j = 0; j < nj+1; ++j) {
		  for(int i = 0; i < ni; ++i) {
			 if(j < 1 || j >= nj || (nodal_solid_phi(i+1,j,k) + nodal_solid_phi(i,j,k)
										+nodal_solid_phi(i+1,j,k+1) + nodal_solid_phi(i,j,k+1))/4.0 <= 0)
				v_state(i,j,k) = SOLID;
			 else
				v_state(i,j,k) = FLUID;
		  }
	   }
   }

   // Determine the solid boundaries for the z-faces. 
   #pragma omp parallel for
   for(int k = 0; k < nk+1; ++k) {
	   for(int j = 0; j < nj; ++j) {
		  for(int i = 0; i < ni; ++i) {
			 if(k < 1 || k >= nk || (nodal_solid_phi(i,j+1,k) + nodal_solid_phi(i,j,k)
										+nodal_solid_phi(i+1,j+1,k) + nodal_solid_phi(i+1,j,k))/4.0 <= 0)
				w_state(i,j,k) = SOLID;
			 else
				w_state(i,j,k) = FLUID;
		  }
	   }
   }

   //printf("Building matrix\n");
   int elts = (ni+1)*nj*nk + ni*(nj+1)*nk + ni*nj*(nk+1);
   if(vrhs.size() != elts) {
      vrhs.resize(elts);
      velocities.resize(elts);
      vmatrix.resize(elts);
   }
   vmatrix.zero();
   vrhs.assign(vrhs.size(), 0);
   velocities.assign(velocities.size(), 0);

   float factor = dt/sqr(dx);
   #pragma omp parallel for
   for(int k = 1; k < nk-1; ++k) for(int j = 1; j < nj-1; ++j) for(int i = 1; i < ni-1; ++i) {
      if(u_state(i,j,k) == FLUID ) {
         int index = u_ind(i,j,k);  
         
         vrhs[index] = u_vol(i,j,k) * u(i,j,k);
         vmatrix.set_element(index,index,u_vol(i,j,k));
         
         //uxx terms
         float visc_right = viscosity(i,j,k);
         float visc_left = viscosity(i-1,j,k);
         float vol_right = c_vol(i,j,k);
         float vol_left = c_vol(i-1,j,k);

        //u_x_right
         vmatrix.add_to_element(index,index, 2*factor*visc_right*vol_right);
         if(u_state(i+1,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i+1, j,k), -2*factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            vrhs[index] -= -2*factor*visc_right*vol_right*u_obj;

         //u_x_left
         vmatrix.add_to_element(index,index, 2*factor*visc_left*vol_left);
         if(u_state(i-1,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i-1,j,k), -2*factor*visc_left*vol_left);
         else if(u_state(i-1,j,k) == SOLID)
            vrhs[index] -= -2*factor*visc_left*vol_left*u_obj;
         
         //uyy terms
         float visc_top = 0.25f*(viscosity(i-1,j+1,k) + viscosity(i-1,j,k) + viscosity(i,j+1,k) + viscosity(i,j,k));
         float visc_bottom = 0.25f*(viscosity(i-1,j,k) + viscosity(i-1,j-1,k) + viscosity(i,j,k) + viscosity(i,j-1,k));
         float vol_top = 0.5f * (n_vol(i,j+1,k) + n_vol(i, j+1, k+1));
         float vol_bottom = 0.5f * (n_vol(i,j,k) + n_vol(i, j, k+1));

         //u_y_top
         vmatrix.add_to_element(index,index, +factor*visc_top*vol_top);
         if(u_state(i,j+1,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i, j+1, k), -factor*visc_top*vol_top);
         else if(u_state(i,j+1,k) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_top*vol_top;
      
         //u_y_bottom
         vmatrix.add_to_element(index,index, +factor*visc_bottom*vol_bottom);
         if(u_state(i,j-1,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j-1,k), -factor*visc_bottom*vol_bottom);
         else if(u_state(i,j-1,k) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_bottom*vol_bottom;

		 //uzz terms
		 float visc_far = 0.25f*(viscosity(i-1,j,k+1) + viscosity(i-1,j,k) + viscosity(i,j,k+1) + viscosity(i,j,k));
		 float visc_near = 0.25f*(viscosity(i-1,j,k) + viscosity(i-1,j,k-1) + viscosity(i,j,k) + viscosity(i,j,k-1));
		 float vol_far = 0.5f * (n_vol(i,j,k+1) + n_vol(i, j+1, k+1));
		 float vol_near = 0.5f * (n_vol(i,j,k) + n_vol(i, j+1, k));

		  //u_z_far
         vmatrix.add_to_element(index,index, +factor*visc_far*vol_far);
         if(u_state(i,j,k+1) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j,k+1), -factor*visc_far*vol_far);
         else if(u_state(i,j,k+1) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_far*vol_far;
      
         //u_z_near
         vmatrix.add_to_element(index,index, +factor*visc_near*vol_near);
         if(u_state(i,j,k-1) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j,k-1), -factor*visc_near*vol_near);
         else if(u_state(i,j,k-1) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_near*vol_near;
      
         //vxy terms
         //v_x_top
         if(v_state(i,j+1,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j+1,k), -factor*visc_top*vol_top);
         else if(v_state(i,j+1,k) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_top*vol_top;
         
         if(v_state(i-1,j+1,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i-1,j+1,k), factor*visc_top*vol_top);
         else if(v_state(i-1,j+1,k) == SOLID)
            vrhs[index] -= v_obj*factor*visc_top*vol_top;
     
         //v_x_bottom
         if(v_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j,k), +factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k) == SOLID)
            vrhs[index] -= v_obj*factor*visc_bottom*vol_bottom;
         
         if(v_state(i-1,j,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i-1,j,k), -factor*visc_bottom*vol_bottom);
         else if(v_state(i-1,j,k) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_bottom*vol_bottom;
      
		 //v_x_far
         if(w_state(i,j,k+1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j,k+1), -factor*visc_far*vol_far);
         else if(w_state(i,j,k+1) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_far*vol_far;
         
         if(w_state(i-1,j,k+1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i-1,j,k+1), factor*visc_far*vol_far);
         else if(w_state(i-1,j,k+1) == SOLID)
            vrhs[index] -= w_obj*factor*visc_far*vol_far;
     
         //v_x_near
         if(w_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j,k), +factor*visc_near*vol_near);
         else if(w_state(i,j,k) == SOLID)
            vrhs[index] -= w_obj*factor*visc_near*vol_near;
         
         if(w_state(i-1,j,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i-1,j,k), -factor*visc_near*vol_near);
         else if(w_state(i-1,j,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_near*vol_near;

      }
   }

   #pragma omp parallel for
   for(int k = 1; k < nk-1; ++k) for(int j = 1; j < nj-1; ++j) for(int i = 1; i < ni-1; ++i) {
      if(v_state(i,j,k) == FLUID ) {
         int index = v_ind(i,j,k);  
         
         vrhs[index] = v_vol(i,j,k) * v(i,j,k);
         vmatrix.set_element(index,index,v_vol(i,j,k));
         
         //vyy terms
         float visc_top = viscosity(i,j,k);
         float visc_bottom = viscosity(i,j-1,k);
         float vol_top = c_vol(i,j,k);
         float vol_bottom = c_vol(i,j-1,k);

        //v_y_top
         vmatrix.add_to_element(index,index, 2*factor*visc_top*vol_top);
         if(v_state(i,j+1,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i, j+1,k), -2*factor*visc_top*vol_top);
         else if(v_state(i,j+1,k) == SOLID)
            vrhs[index] -= -2*factor*visc_top*vol_top*v_obj;

         //v_y_bottom
         vmatrix.add_to_element(index,index, 2*factor*visc_bottom*vol_bottom);
         if(v_state(i,j-1,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j-1,k), -2*factor*visc_bottom*vol_bottom);
         else if(v_state(i,j-1,k) == SOLID)
            vrhs[index] -= -2*factor*visc_bottom*vol_bottom*v_obj;
         
         //vxx terms
         float visc_right = 0.25f*(viscosity(i+1,j-1,k) + viscosity(i,j-1,k) + viscosity(i+1,j,k) + viscosity(i,j,k));
         float visc_left = 0.25f*(viscosity(i,j-1,k) + viscosity(i-1,j-1,k) + viscosity(i,j,k) + viscosity(i-1,j,k));
         float vol_right = 0.5f * (n_vol(i+1,j,k) + n_vol(i+1, j, k+1));
         float vol_left = 0.5f * (n_vol(i,j,k) + n_vol(i, j, k+1));

         //v_x_right
         vmatrix.add_to_element(index,index, +factor*visc_right*vol_right);
         if(v_state(i+1,j,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i+1, j, k), -factor*visc_right*vol_right);
         else if(v_state(i+1,j,k) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_right*vol_right;
      
         //v_x_left
         vmatrix.add_to_element(index,index, +factor*visc_left*vol_left);
         if(v_state(i-1,j,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i-1,j,k), -factor*visc_left*vol_left);
         else if(v_state(i-1,j,k) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_left*vol_left;

		 //vzz terms
		 float visc_far = 0.25f*(viscosity(i,j-1,k) + viscosity(i,j-1,k+1) + viscosity(i,j,k) + viscosity(i,j,k+1));
		 float visc_near = 0.25f*(viscosity(i,j-1,k) + viscosity(i,j-1,k-1) + viscosity(i,j,k) + viscosity(i,j,k-1));
		 float vol_far = 0.5f * (n_vol(i,j,k+1) + n_vol(i+1, j, k+1));
		 float vol_near = 0.5f * (n_vol(i,j,k) + n_vol(i+1, j, k));

		  //u_z_far
         vmatrix.add_to_element(index,index, +factor*visc_far*vol_far);
         if(v_state(i,j,k+1) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j,k+1), -factor*visc_far*vol_far);
         else if(v_state(i,j,k+1) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_far*vol_far;
      
         //u_z_near
         vmatrix.add_to_element(index,index, +factor*visc_near*vol_near);
         if(v_state(i,j,k-1) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j,k-1), -factor*visc_near*vol_near);
         else if(v_state(i,j,k-1) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_near*vol_near;
      
         //vxy terms
         //right
         if(u_state(i+1,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i+1,j,k), -factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_right*vol_right;
         
         if(u_state(i+1,j-1,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i+1,j-1,k), factor*visc_right*vol_right);
         else if(u_state(i+1,j-1,k) == SOLID)
            vrhs[index] -= u_obj*factor*visc_right*vol_right;
     
         //left
         if(u_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j,k), +factor*visc_left*vol_left);
         else if(u_state(i,j,k) == SOLID)
            vrhs[index] -= u_obj*factor*visc_left*vol_left;
         
         if(u_state(i,j-1,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j-1,k), -factor*visc_left*vol_left);
         else if(u_state(i,j-1,k) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_left*vol_left;
      
		 //v_far
         if(w_state(i,j,k+1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j,k+1), -factor*visc_far*vol_far);
         else if(w_state(i,j,k+1) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_far*vol_far;
         
         if(w_state(i,j-1,k+1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j-1,k+1), factor*visc_far*vol_far);
         else if(w_state(i,j-1,k+1) == SOLID)
            vrhs[index] -= w_obj*factor*visc_far*vol_far;
     
         //near
         if(w_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j,k), +factor*visc_near*vol_near);
         else if(w_state(i,j,k) == SOLID)
            vrhs[index] -= w_obj*factor*visc_near*vol_near;
         
         if(w_state(i,j-1,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j-1,k), -factor*visc_near*vol_near);
         else if(w_state(i,j-1,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_near*vol_near;

      }
   }

   #pragma omp parallel for
   for(int k = 1; k < nk-1; ++k) for(int j = 1; j < nj-1; ++j) for(int i = 1; i < ni-1; ++i) {
      if(w_state(i,j,k) == FLUID ) {
         int index = w_ind(i,j,k);  
         
         vrhs[index] = w_vol(i,j,k) * w(i,j,k);
         vmatrix.set_element(index,index,w_vol(i,j,k));
         
         //vyy terms
         float visc_far = viscosity(i,j,k);
         float visc_near = viscosity(i,j,k-1);
         float vol_far = c_vol(i,j,k);
         float vol_near = c_vol(i,j,k-1);

        //v_y_top
         vmatrix.add_to_element(index,index, 2*factor*visc_far*vol_far);
         if(w_state(i,j,k+1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i, j,k+1), -2*factor*visc_far*vol_far);
         else if(w_state(i,j,k+1) == SOLID)
            vrhs[index] -= -2*factor*visc_far*vol_far*w_obj;

         //v_y_bottom
         vmatrix.add_to_element(index,index, 2*factor*visc_near*vol_near);
         if(w_state(i,j,k-1) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j,k-1), -2*factor*visc_near*vol_near);
         else if(w_state(i,j,k-1) == SOLID)
            vrhs[index] -= -2*factor*visc_near*vol_near*w_obj;
         
         //vxx terms
         float visc_right = 0.25f*(viscosity(i,j,k-1) + viscosity(i+1,j,k-1) + viscosity(i,j,k) + viscosity(i+1,j,k));
         float visc_left = 0.25f*(viscosity(i,j,k-1) + viscosity(i-1,j,k-1) + viscosity(i,j,k) + viscosity(i-1,j,k));
         float vol_right = 0.5f * (n_vol(i+1,j,k) + n_vol(i+1, j+1, k));
         float vol_left = 0.5f * (n_vol(i,j,k) + n_vol(i, j+1, k));

         //v_x_right
         vmatrix.add_to_element(index,index, +factor*visc_right*vol_right);
         if(w_state(i+1,j,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i+1, j, k), -factor*visc_right*vol_right);
         else if(w_state(i+1,j,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_right*vol_right;
      
         //v_x_left
         vmatrix.add_to_element(index,index, +factor*visc_left*vol_left);
         if(w_state(i-1,j,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i-1,j,k), -factor*visc_left*vol_left);
         else if(w_state(i-1,j,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_left*vol_left;

		 //vzz terms
		 float visc_top = 0.25f*(viscosity(i,j,k-1) + viscosity(i,j+1,k-1) + viscosity(i,j,k) + viscosity(i,j+1,k));
		 float visc_bottom = 0.25f*(viscosity(i,j,k-1) + viscosity(i,j-1,k-1) + viscosity(i,j,k) + viscosity(i,j-1,k));
		 float vol_top = 0.5f * (n_vol(i,j+1,k) + n_vol(i+1, j+1, k));
		 float vol_bottom = 0.5f * (n_vol(i,j,k) + n_vol(i+1, j, k));

		  //u_z_far
         vmatrix.add_to_element(index,index, +factor*visc_top*vol_top);
         if(w_state(i,j+1,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j+1,k), -factor*visc_top*vol_top);
         else if(w_state(i,j+1,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_top*vol_top;
      
         //u_z_near
         vmatrix.add_to_element(index,index, +factor*visc_bottom*vol_bottom);
         if(w_state(i,j-1,k) == FLUID)
            vmatrix.add_to_element(index,w_ind(i,j-1,k), -factor*visc_bottom*vol_bottom);
         else if(w_state(i,j-1,k) == SOLID)
            vrhs[index] -= -w_obj*factor*visc_bottom*vol_bottom;
      
         //vxy terms
         //top
         if(v_state(i,j+1,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j+1,k), -factor*visc_top*vol_top);
         else if(v_state(i,j+1,k) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_top*vol_top;
         
         if(v_state(i,j+1,k-1) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j+1,k-1), factor*visc_top*vol_top);
         else if(v_state(i,j+1,k-1) == SOLID)
            vrhs[index] -= v_obj*factor*visc_top*vol_top;
     
         //bottom
         if(v_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j,k), +factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k) == SOLID)
            vrhs[index] -= v_obj*factor*visc_bottom*vol_bottom;
         
         if(v_state(i,j,k-1) == FLUID)
            vmatrix.add_to_element(index,v_ind(i,j,k-1), -factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k-1) == SOLID)
            vrhs[index] -= -v_obj*factor*visc_bottom*vol_bottom;
      
		 //right
         if(u_state(i+1,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i+1,j,k), -factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_right*vol_right;
			
         if(u_state(i+1,j,k-1) == FLUID)
            vmatrix.add_to_element(index,u_ind(i+1,j,k-1), factor*visc_right*vol_right);
         else if(u_state(i+1,j,k-1) == SOLID)
            vrhs[index] -= u_obj*factor*visc_right*vol_right;
     
         //left
         if(u_state(i,j,k) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j,k), +factor*visc_left*vol_left);
         else if(u_state(i,j,k) == SOLID)
            vrhs[index] -= u_obj*factor*visc_left*vol_left;
         
         if(u_state(i,j,k-1) == FLUID)
            vmatrix.add_to_element(index,u_ind(i,j,k-1), -factor*visc_left*vol_left);
         else if(u_state(i,j,k-1) == SOLID)
            vrhs[index] -= -u_obj*factor*visc_left*vol_left;

      }
   }

   double res_out;
   int iter_out;

   //std::cout << "about to solve" << std::endl;
   solver.set_solver_parameters(1e-18, 1000);
   solver.solve(vmatrix, vrhs, velocities, res_out, iter_out);
   //std::cout << "done solve" << std::endl;

   // Update the final velocities like a boss. 
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) {
    for(int j = 0; j < nj; ++j) {
      for(int i = 0; i < ni+1; ++i) {
         if(u_state(i,j,k) == FLUID)
            u(i,j,k) = (float)velocities[u_ind(i,j,k)];
         else if(u_state(i,j,k) == SOLID) 
            u(i,j,k) = u_obj;
	  }
	}
   }
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) {
    for(int j = 0; j < nj+1; ++j) {
      for(int i = 0; i < ni; ++i) {
         if(v_state(i,j,k) == FLUID)
            v(i,j,k) = (float)velocities[v_ind(i,j,k)];
         else if(v_state(i,j,k) == SOLID) 
            v(i,j,k) = v_obj;
	  }
	}
   }
   #pragma omp parallel for
   for(int k = 0; k < nk+1; ++k) {
    for(int j = 0; j < nj; ++j) {
      for(int i = 0; i < ni; ++i) {
         if(w_state(i,j,k) == FLUID)
            w(i,j,k) = (float)velocities[w_ind(i,j,k)];
         else if(w_state(i,j,k) == SOLID) 
            w(i,j,k) = w_obj;
	  }
	}
   }
}

void FluidSim::compute_volume_fractions(const Array3f& levelset, Array3f& fractions, glm::vec3 fraction_origin, int subdivision) {
	//Assumes levelset and fractions have the same dx
   float sub_dx = 1.0 / subdivision;
   int sample_max = subdivision*subdivision*subdivision;

   #pragma omp parallel for
   for(int k = 0; k < fractions.nk; ++k) {
	   for(int j = 0; j < fractions.nj; ++j) {
		  for(int i = 0; i < fractions.ni; ++i) {
			 float start_x = fraction_origin[0] + (float)i;
			 float start_y = fraction_origin[1] + (float)j;
			 float start_z = fraction_origin[2] + (float)k;
			 int incount = 0;

			 for(int sub_k = 0; sub_k < subdivision; ++sub_k) {
				 for(int sub_j = 0; sub_j < subdivision; ++sub_j) {
					for(int sub_i = 0; sub_i < subdivision; ++sub_i) {
					   float x_pos = start_x + (sub_i+0.5)*sub_dx;
					   float y_pos = start_y + (sub_j+0.5)*sub_dx;
					   float z_pos = start_z + (sub_k+0.5)*sub_dx;
					   float phi_val = interpolate_value<glm::vec3>(glm::vec3(x_pos,y_pos,z_pos), levelset);
					   if(phi_val < 0) 
						  ++incount;
					}
				 }
			 }
			 fractions(i,j,k) = (float)incount / (float)sample_max;
		  }
	   }
   }
}

void FluidSim::compute_viscosity_weights(float dt) {
   compute_volume_fractions(liquid_phi, c_vol, glm::vec3(-0.5,-0.5,-0.5), 2);
   compute_volume_fractions(liquid_phi, n_vol, glm::vec3(-1, -1,-1), 2);
   compute_volume_fractions(liquid_phi, u_vol, glm::vec3(-1,-0.5,-0.5), 2);
   compute_volume_fractions(liquid_phi, v_vol, glm::vec3(-0.5,-1,-0.5), 2);
   compute_volume_fractions(liquid_phi, w_vol, glm::vec3(-0.5,-0.5,-1), 2);
}

void FluidSim::project(float dt) {

   //Estimate the liquid signed distance
   compute_phi();
   
   //Compute finite-volume type face area weight for each velocity sample.
   compute_weights();

   //Set up and solve the variational pressure solve.
   solve_pressure(dt);
   
}


//Apply RK2 to advect a point in the domain.
glm::vec3 FluidSim::trace_rk2(const glm::vec3& position, float dt) {
   glm::vec3 input = position;
   glm::vec3 velocity = get_velocity(input);
   velocity = get_velocity(input + 0.5f*dt*velocity);
   input += dt*velocity;
   return input;
}

//Interpolate velocity from the MAC grid.
glm::vec3 FluidSim::get_velocity(const glm::vec3& position) {

   //Interpolate the velocity from the u and v grids
    float u_value = interpolate_value<glm::vec3>(position / dx - glm::vec3(0, 0.5f, 0.5f), u);
    float v_value = interpolate_value<glm::vec3>(position / dx - glm::vec3(0.5f, 0, 0.5f), v);
    float w_value = interpolate_value<glm::vec3>(position / dx - glm::vec3(0.5f, 0.5f, 0), w);

    return glm::vec3(u_value, v_value, w_value);
}



//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim::compute_weights() {

   //Compute face area fractions (using marching squares cases).
	#pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      u_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,j,  k),
                                             nodal_solid_phi(i,j+1,k),
                                             nodal_solid_phi(i,j,  k+1),
                                             nodal_solid_phi(i,j+1,k+1));
      u_weights(i,j,k) = clamp(u_weights(i,j,k),0.0f,1.0f);
   }
   #pragma omp parallel for
   for(int k = 0; k < nk; ++k) for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      v_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,k),
                                             nodal_solid_phi(i,  j,k+1),
                                             nodal_solid_phi(i+1,j,k),
                                             nodal_solid_phi(i+1,j,k+1));
      v_weights(i,j,k) = clamp(v_weights(i,j,k),0.0f,1.0f);
   }
   #pragma omp parallel for
   for(int k = 0; k < nk+1; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      w_weights(i,j,k) = 1 - fraction_inside(nodal_solid_phi(i,  j,  k),
                                             nodal_solid_phi(i,  j+1,k),
                                             nodal_solid_phi(i+1,j,  k),
                                             nodal_solid_phi(i+1,j+1,k));
      w_weights(i,j,k) = clamp(w_weights(i,j,k),0.0f,1.0f);
   }


}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim::solve_pressure(float dt) {


   int ni = v.ni;
   int nj = u.nj;
   int nk = u.nk;

   int system_size = ni*nj*nk;
   if(rhs.size() != system_size) {
      rhs.resize(system_size);
      pressure.resize(system_size);
      matrix.resize(system_size);
   }
   
   matrix.zero();
   rhs.assign(rhs.size(), 0);
   pressure.assign(pressure.size(), 0);

   #pragma omp parallel for
   //Build the linear system for pressure
   for(int k = 1; k < nk-1; ++k) {
      for(int j = 1; j < nj-1; ++j) {
         for(int i = 1; i < ni-1; ++i) {
            int index = i + ni*j + ni*nj*k;

            rhs[index] = 0;
            pressure[index] = 0;
            float centre_phi = liquid_phi(i,j,k);
            if(centre_phi < 0) {

               //right neighbour
               float term = u_weights(i+1,j,k) * dt / sqr(dx);
               float right_phi = liquid_phi(i+1,j,k);
               if(right_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + 1, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, right_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= u_weights(i+1,j,k)*u(i+1,j,k) / dx;

               //left neighbour
               term = u_weights(i,j,k) * dt / sqr(dx);
               float left_phi = liquid_phi(i-1,j,k);
               if(left_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - 1, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, left_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += u_weights(i,j,k)*u(i,j,k) / dx;

               //top neighbour
               term = v_weights(i,j+1,k) * dt / sqr(dx);
               float top_phi = liquid_phi(i,j+1,k);
               if(top_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, top_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= v_weights(i,j+1,k)*v(i,j+1,k) / dx;

               //bottom neighbour
               term = v_weights(i,j,k) * dt / sqr(dx);
               float bot_phi = liquid_phi(i,j-1,k);
               if(bot_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, bot_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += v_weights(i,j,k)*v(i,j,k) / dx;


               //far neighbour
               term = w_weights(i,j,k+1) * dt / sqr(dx);
               float far_phi = liquid_phi(i,j,k+1);
               if(far_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, far_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= w_weights(i,j,k+1)*w(i,j,k+1) / dx;

               //near neighbour
               term = w_weights(i,j,k) * dt / sqr(dx);
               float near_phi = liquid_phi(i,j,k-1);
               if(near_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, near_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += w_weights(i,j,k)*w(i,j,k) / dx;

               /*
               //far neighbour
               term = w_weights(i,j,k+1) * dt / sqr(dx);
               float far_phi = liquid_phi(i,j,k+1);
               if(far_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index + ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, far_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] -= w_weights(i,j,k+1)*w(i,j,k+1) / dx;

               //near neighbour
               term = w_weights(i,j,k) * dt / sqr(dx);
               float near_phi = liquid_phi(i,j,k-1);
               if(near_phi < 0) {
                  matrix.add_to_element(index, index, term);
                  matrix.add_to_element(index, index - ni*nj, -term);
               }
               else {
                  float theta = fraction_inside(centre_phi, near_phi);
                  if(theta < 0.01f) theta = 0.01f;
                  matrix.add_to_element(index, index, term/theta);
               }
               rhs[index] += w_weights(i,j,k)*w(i,j,k) / dx;   
               */

            }
         }
      }
   }

   //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

   double tolerance;
   int iterations;
   solver.set_solver_parameters(1e-18, 1000);
   bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
   //printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
   if(!success) {
      printf("WARNING: Pressure solve failed!************************************************\n");
   }

   //Apply the velocity update
   u_valid.assign(0);
   #pragma omp parallel for
   for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(u_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i-1,j,k) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
            theta = fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         u(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-1]) / dx / theta; 
         u_valid(i,j,k) = 1;
      }
   }
   
   v_valid.assign(0);
   #pragma omp parallel for
   for(int k = 0; k < v.nk; ++k) for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(v_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j-1,k) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
            theta = fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         v(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni]) / dx / theta; 
         v_valid(i,j,k) = 1;
      }
   }

   w_valid.assign(0);
   #pragma omp parallel for
   for(int k = 1; k < w.nk-1; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      int index = i + j*ni + k*ni*nj;
      if(w_weights(i,j,k) > 0 && (liquid_phi(i,j,k) < 0 || liquid_phi(i,j,k-1) < 0)) {
         float theta = 1;
         if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
            theta = fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
         if(theta < 0.01f) theta = 0.01f;
         w(i,j,k) -= dt  * (float)(pressure[index] - pressure[index-ni*nj]) / dx / theta; 
         w_valid(i,j,k) = 1;
      }
   }
 

   for(unsigned int i = 0; i < u_valid.a.size(); ++i)
      if(u_valid.a[i] == 0)
         u.a[i] = 0;
   for(unsigned int i = 0; i < v_valid.a.size(); ++i)
      if(v_valid.a[i] == 0)
         v.a[i] = 0;
   for(unsigned int i = 0; i < w_valid.a.size(); ++i)
      if(w_valid.a[i] == 0)
         w.a[i] = 0;
}


//Apply several iterations of a very simple propagation of valid velocity data in all directions
void extrapolate(Array3f& grid, Array3c& valid) {

   Array3f temp_grid = grid;
   Array3c old_valid(valid.ni,valid.nj,valid.nk);
   #pragma omp parallel for
   for(int layers = 0; layers < 10; ++layers) {
      old_valid = valid;
      for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) {
         float sum = 0;
         int count = 0;

         if(!old_valid(i,j,k)) {

            if(old_valid(i+1,j,k)) {
               sum += grid(i+1,j,k);
               ++count;
            }
            if(old_valid(i-1,j,k)) {
               sum += grid(i-1,j,k);
               ++count;
            }
            if(old_valid(i,j+1,k)) {
               sum += grid(i,j+1,k);
               ++count;
            }
            if(old_valid(i,j-1,k)) {
               sum += grid(i,j-1,k);
               ++count;
            }
            if(old_valid(i,j,k+1)) {
               sum += grid(i,j,k+1);
               ++count;
            }
            if(old_valid(i,j,k-1)) {
               sum += grid(i,j,k-1);
               ++count;
            }

            //If any of neighbour cells were valid, 
            //assign the cell their average value and tag it as valid
            if(count > 0) {
               temp_grid(i,j,k) = sum /(float)count;
               valid(i,j,k) = 1;
            }

         }
      }
      grid = temp_grid;

   }

}

void FluidSim::draw() {

//	glEnable(GL_LIGHTING);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//   set_lights_and_material(0); 
    
    //-------------------------------------------------------------
    //-------------------DRAW FLUID PARTICLES----------------------
    //-------------------------------------------------------------
    
    
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
   GLUquadric* particle_sphere;
   particle_sphere = gluNewQuadric();
   gluQuadricDrawStyle(particle_sphere, GLU_FILL );

   for(unsigned int p = 0; p < particles.size(); ++p) {
	   
      glPushMatrix();
	  glm::vec3 pos = particles[p]->position;
      glTranslatef(pos[0]-(.5*grid_width), pos[1]-(.5*grid_width), pos[2]-(.5*grid_width));
	  gluQuadricNormals(particle_sphere, GLU_SMOOTH);
	  if(transparentRender){
		   glColor4f(particles[p]->color[0], particles[p]->color[1], particles[p]->color[2], 0.2);
	  }else{
		   glColor4f(particles[p]->color[0], particles[p]->color[1], particles[p]->color[2], 1.0);
	  }
	 
      gluSphere(particle_sphere, particle_radius, 20, 20);
      glPopMatrix();   
   }

    //-------------------------------------------------------------
    //-----------------DRAW GRID BOUNDING BOX----------------------
    //-------------------------------------------------------------
    
    //calculate boundary points
    float p = .5 * grid_width;
    
  // glDisable(GL_LIGHTING);
	glColor3f(0,0,0);
	glBegin(GL_LINES);
	glVertex3f(-p,-p,-p);
	glVertex3f(-p,-p,p);
    
	glVertex3f(-p,-p,-p);
	glVertex3f(-p,p,-p);
    
	glVertex3f(-p,-p,-p);
	glVertex3f(p,-p,-p);
    
	glVertex3f(-p,p,-p);
	glVertex3f(p,p,-p);
    
	glVertex3f(p,p,-p);
	glVertex3f(p,-p,-p);
    
	glVertex3f(p,-p,-p);
	glVertex3f(p,-p,p);
    
	glVertex3f(-p,p,-p);
	glVertex3f(-p,p,p);
    
	glVertex3f(p,p,-p);
	glVertex3f(p,p,p);
    
	glVertex3f(-p,p,p);
	glVertex3f(p,p,p);
    
	glVertex3f(p,-p,p);
	glVertex3f(p,p,p);
    
	glVertex3f(-p,-p,p);
	glVertex3f(p,-p,p);
    
	glVertex3f(-p,-p,p);
	glVertex3f(-p,p,p);
	glEnd();
   
   //Draw wireframe sphere geometry (specific to this scene).
   glColor3f(0,0,0);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
   GLUquadric* sphere;
   sphere = gluNewQuadric();
   gluQuadricDrawStyle(sphere, GLU_LINE );
   glPushMatrix();
  // glTranslatef(0.0f, 0.0f,0.0f);
	//gluSphere(sphere, 0.30, 20, 20);
   //glTranslatef(0.55f-.5, 0.8f-.5, 0.4f-.5);
   //gluSphere(sphere, 0.23f, 20, 20);
   glPopMatrix();
    

   
   /*for(int x = 0; x<mesh->getObjBounds().x; x++) {
	   for(int y = 0; y<mesh->getObjBounds().y; y++) {
		   for(int z = 0; z<mesh->getObjBounds().z; z++) {
			   
			 if(mesh->getSignedDistanceAtCellObjSpace(x,y,z)<.0 && mesh->getSignedDistanceAtCellObjSpace(x,y,z)>-100){
				   glColor3f(0,0,0);
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
				GLUquadric* sphere;
				sphere = gluNewQuadric();
				gluQuadricDrawStyle(sphere, GLU_LINE );
				glPushMatrix();
				float xt = ((x)-(.5*mesh->getObjBounds().x))/(mesh->getScaleFactor());
				float yt = ((y)-(.5*mesh->getObjBounds().y))/(mesh->getScaleFactor());
				float zt = ((z)-(.5*mesh->getObjBounds().z))/(mesh->getScaleFactor());
				//std::cout << xt<< " " << yt<< " " << zt <<std::endl;
				glTranslatef(xt, yt, zt);
				gluSphere(sphere, 0.005, 2, 2);
				glPopMatrix();
			 }
		   }
	   }
   }*/

    /* for(int x = 0; x<mesh->getObjBounds().x; x++) {
	   for(int y = 0; y<mesh->getObjBounds().y; y++) {
		   for(int z = 0; z<mesh->getObjBounds().z; z++) {
			   	float xt = ((x)-(.5*mesh->getObjBounds().x))/(mesh->getScaleFactor());
				float yt = ((y)-(.5*mesh->getObjBounds().y))/(mesh->getScaleFactor());
				float zt = ((z)-(.5*mesh->getObjBounds().z))/(mesh->getScaleFactor());
			 if(mesh->getSignedDistanceAtCellWorldSpace(xt,yt,zt)<.03 && mesh->getSignedDistanceAtCellWorldSpace(xt,yt,zt)>-.03){
				   glColor3f(1,0,0);
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
				GLUquadric* sphere;
				sphere = gluNewQuadric();
				gluQuadricDrawStyle(sphere, GLU_LINE );
				glPushMatrix();

				glTranslatef(xt, yt, zt);
				gluSphere(sphere, 0.005, 2, 2);
				glPopMatrix();
			 }
		   }
	   }
   }*/

/*glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
for(int i=0; i<mesh->faceList.size(); i++){
	glm::vec3 face = mesh->faceList[i];
	glBegin(GL_LINES);

	glm::vec3 point0 = mesh->vertList[face[0]]-mesh->meshCenter ;
	glm::vec3 point1 = mesh->vertList[face[1]]-mesh->meshCenter ;
	glm::vec3 point2 = mesh->vertList[face[2]]-mesh->meshCenter ;

	glVertex3f(point0[0],point0[1],point0[2]);
	glVertex3f(point1[0],point1[1],point1[2]);

	glVertex3f(point1[0],point1[1],point1[2]);
	glVertex3f(point2[0],point2[1],point2[2]);

	glVertex3f(point2[0],point2[1],point2[2]);
	glVertex3f(point0[0],point0[1],point0[2]);
	glEnd();
}

glPopMatrix();*/

   //run marching cubes for mesh reconstruction
   //MarchingCubes(liquid_phi, dx, frameNum, outputOBJ);

   frameNum++;

   if (mRecordEnabled) grabScreen();
}

void FluidSim::set_lights_and_material(int object)
{
   glEnable(GL_LIGHTING);
   GLfloat global_ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
   glShadeModel(GL_SMOOTH);

   //Light #1
   GLfloat color[4] = {0.0f, 0.0f, 1.0f, 1.0f};
   GLfloat position[3] = {1.0f, 1.0f, 1.0f};
  // glLightfv(GL_LIGHT0, GL_SPECULAR, color);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
   glLightfv(GL_LIGHT0, GL_POSITION, position);

   //Light #2
   GLfloat color2[4] = {0.0f, 0.0f, 1.0f, 1.0f};
   GLfloat position2[3] = {-1.0f, -1.0f, 1.0f};
  // glLightfv(GL_LIGHT1, GL_SPECULAR, color2);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, color2);
   glLightfv(GL_LIGHT1, GL_POSITION, position2);

   GLfloat obj_color[4] = {.2, .3, .7};
  // glMaterialfv (GL_FRONT, GL_AMBIENT, obj_color);
   glMaterialfv (GL_FRONT, GL_DIFFUSE, obj_color);

  // GLfloat specular[4] = {.4, .2, .8};
  // glMaterialf (GL_FRONT, GL_SHININESS, 32);
  // glMaterialfv (GL_FRONT, GL_SPECULAR, specular);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);
}
