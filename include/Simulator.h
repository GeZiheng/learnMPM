#pragma once
#include <vector>
#include "Particle.h"
#include "Grid.h"
#include "ParticleSystem.h"

class Simulator
{
private:
	ParticleSystem		pts_sys;			// particles system
	ParticleSystem		pts_sys_pre;		// particles system prediction
	double				t;					// time

private:
	void predictor();						// predictor step
	void corrector();						// corrector step
	void writeData(int frame_num);			// write particles data into file
	void writeDataa(int frame_num);			// write particles data into file
	void renderScene();						// render scene on the screen
	void densityRecons();					// density reconstruction

public:
	Simulator();
	~Simulator();
	void init();							// initialize parameters and sample particles
	void oneTimeStep();						// execute a whole cycle in one time step
	double getTime();						// return current time
};