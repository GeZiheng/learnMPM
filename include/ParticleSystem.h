#pragma once
#include "Particle.h"
#include "Grid.h"

class ParticleSystem
{
public:
	vector<Particle>	particles;						// particle cloud
	int					n_real;							// number of real particles
	Grid				grid;							// grid for neighbor search

public:
	ParticleSystem();
	~ParticleSystem();
	ParticleSystem(int n, int d);
	void addParticleToGrid(Particle particle, int id);	// add one particle to the grid
	void countParticlesOnGrid();						// count how many particles are inside each grid node
	vector<int> findNeighbors(int id, double radius);	// find neighbors of one particle
	void addGhostParticles();							// add ghost particles
	void deleteGhostParticles();						// delete ghost particles
	void handleCollisions();							// correct particle positions for collisions
};