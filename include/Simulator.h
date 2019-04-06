#pragma once
#include <vector>
#include "Particle.h"
#include "Grid.h"

class Simulator
{
private:
	int					pts_num;			// number of particles
	vector<Particle>	pts_cloud;			// particle cloud
	Grid*				grid;				// grid
	double				t;					// time

public:
	double				dx;					// grid dx width
	double				dt;					// time step length
	double				rho0;				// initial density
	double				mu;					// mu
	double				lambda;				// lambda
	double				gravity;			// gravity
	Vector2d			rec_center;			// center position
	double				rec_size;			// rectangle size
	double				end_t;				// end time
	int					frame_num;			// number of frame
	double				frame_dt;			// frame dt

private:
	void P2G();								// particle to grid transfer
	void G2P();								// grid to particle transfer
	void gridV();							// update grid velocity
	void gridF();							// update grid forces
	void ptsDG();							// evolve partivle deformation gradient
	void advection();						// particle advection
	void writeData(int frame_num);			// write particles data into file

public:
	Simulator();
	~Simulator();
	void init();							// initialize parameters and sample particles
	void oneTimeStep();						// execute a whole cycle in one time step
	double getTime();						// return current time
};