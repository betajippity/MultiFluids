#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "array3_utils.h"
#include <vector>
#include "objsdf.h"

struct particle{
	glm::vec3 position;
	glm::vec3 color;
	double viscosity;

	particle(){
		position = glm::vec3(0,0,0);
		color = glm::vec3(0,0,0);
		viscosity = 0.0;
	}
};

class FluidSim {



public:
    
    FluidSim();
    ~FluidSim();
    
    objsdf* mesh;

    void initialize(float width, int ni_, int nj_, int nk_);
    void set_boundary(float (*phi)(const glm::vec3&));
    void set_liquid(float (*phi)(const glm::vec3&), glm::vec3 color, double viscosityCoef);
    void add_particle(const glm::vec3& pos, const glm::vec3& color);

    void advance(float dt);

	void set_lights_and_material(int object);
	void draw();

    //Grid dimensions
    int ni,nj,nk;
    float dx;
    float grid_width;
    
    //Fluid velocity
    Array3f u, v, w;
    Array3f temp_u, temp_v, temp_w;

    //Static geometry representation
    Array3f nodal_solid_phi;
    Array3f u_weights, v_weights, w_weights;
    Array3c u_valid, v_valid, w_valid;


	std::vector<particle*> particles;
    //std::vector<glm::vec3> particles;
	//std::vector<glm::vec3> colors;

    float particle_radius;

    Array3f liquid_phi;

    //Data arrays for extrapolation
    Array3c valid, old_valid;

    //Solver data
    PCGSolver<double> solver;
    SparseMatrixd matrix;
    std::vector<double> rhs;
    std::vector<double> pressure;


	// Viscosity Information
	void apply_viscosity(float dt);
	void compute_viscosity_weights(float dt);
	void compute_volume_fractions(const Array3f& levelset, Array3f& fractions, glm::vec3 fraction_origin, int subdivision);
	void solve_viscosity(float dt);

	int u_ind(int i, int j, int k);
	int v_ind(int i, int j, int k);
	int w_ind(int i, int j, int k);

	//Data for viscosity solve
    Array3f u_vol, v_vol, w_vol, c_vol, n_vol;
    Array3f viscosity;
	SparseMatrixd vmatrix;
    std::vector<double> vrhs;
    std::vector<double> velocities;

    glm::vec3 get_velocity(const glm::vec3& position);

	int getTotalFrames();
	bool isRecording();
	void setRecording(bool on, int width, int height);
	void grabScreen();
	void reset(float width, int ni_, int nj_, int nk_, float (*phi)(const glm::vec3&));
	void setTransparentRender(bool set);
	bool isTransparentRender();
	void setVerbose(bool set);
	bool isVerbose();

	bool verbose;
	bool transparentRender;

	bool outputOBJ;
	int frameNum;

private:

   glm::vec3 trace_rk2(const glm::vec3& position, float dt);

   float cfl();

   void advect_particles(float dt);
   void advect(float dt);
   void add_force(float dt);
   void project(float dt);
   void constrain_velocity();

   //helpers for pressure projection
   void compute_weights();
   void solve_pressure(float dt);
   void compute_phi();

   int mTotalFrameNum;
   bool mRecordEnabled;
   int mFrameNum;
   int recordWidth;
	int recordHeight;
};


#endif