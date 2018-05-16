#include <iostream>
#include "Parameters.h"
#include "ParticleSystem.h"

ParticleSystem::ParticleSystem()
{
}

ParticleSystem::~ParticleSystem()
{
}

ParticleSystem::ParticleSystem(int n, int d)
{
	n_real = 0;
	grid = Grid(n, d);
}

void ParticleSystem::addParticleToGrid(Particle particle, int id)
{
	int grid_i = grid.wid + floor(particle.x.x() / dx);
	int grid_j = grid.wid + floor(particle.x.y() / dx);
	if (grid_i >= 0 && grid_i < grid.num + 2 * grid.wid && grid_j >= 0 && grid_j < grid.num + 2 * grid.wid)
		grid.elem[grid_i][grid_j].push_back(id);
}

void ParticleSystem::countParticlesOnGrid()
{
	int i, j;
	for (i = 0; i < grid.num + 2 * grid.wid; i++)
		for (j = 0; j < grid.num + 2 * grid.wid; j++)
			grid.elem[i][j].clear();
	for (i = 0; i < particles.size(); i++)
		addParticleToGrid(particles[i], i);
};

vector<int> ParticleSystem::findNeighbors(int id, double radius)
{
	int i, j, k;
	vector<int> neighbors;
	Vector2d nei_pos;
	Vector2d pos = particles[id].x;
	double dis;
	for (i = grid.wid + int(floor((pos.x() - radius) / dx)); i <= grid.wid + int(floor((pos.x() + radius) / dx)); i++)
		for (j = grid.wid + int(floor((pos.y() - radius) / dx)); j <= grid.wid + int(floor((pos.y() + radius) / dx)); j++)
			for (k = 0; k < grid.elem[i][j].size(); k++)
			{
				nei_pos = particles[grid.elem[i][j][k]].x;
				dis = (nei_pos - pos).norm();
				if (dis > 0 && dis < radius)
					neighbors.push_back(grid.elem[i][j][k]);
			}
	return neighbors;
}

void ParticleSystem::addGhostParticles()
{
	int i, j, k, l;
	double new_x, new_y, new_vx, new_vy;
	Particle ghost;
	for (i = grid.wid; i < grid.wid + res; i++)
		for (j = 0; j < grid.wid; j++)
		{
			for (k = 0; k < grid.elem[i][grid.wid + j].size(); k++)
			{
				l = grid.elem[i][grid.wid + j][k];
				new_x = particles[l].x.x();
				new_y = - particles[l].x.y();
				new_vx = particles[l].v.x();
				new_vy = - particles[l].v.y();
				ghost = Particle(new_x, new_y, new_vx, new_vy, particles[l].rho, particles[l].p, particles[l].m, Ghost);
				addParticleToGrid(ghost, particles.size());
				particles.push_back(ghost);
			}
			for (k = 0; k < grid.elem[i][grid.wid + res - 1 - j].size(); k++)
			{
				l = grid.elem[i][grid.wid + res - 1 - j][k];
				new_x = particles[l].x.x();
				new_y = 2.0 * dom_size - particles[l].x.y();
				new_vx = particles[l].v.x();
				new_vy = - particles[l].v.y();
				ghost = Particle(new_x, new_y, new_vx, new_vy, particles[l].rho, particles[l].p, particles[l].m, Ghost);
				addParticleToGrid(ghost, particles.size());
				particles.push_back(ghost);
			}
		}
	for (i = 0; i < grid.wid; i++)
		for (j = 0; j < 2 * grid.wid + res; j++)
		{
			for (k = 0; k < grid.elem[grid.wid + i][j].size(); k++)
			{
				l = grid.elem[grid.wid + i][j][k];
				new_x = - particles[l].x.x();
				new_y = particles[l].x.y();
				new_vx = - particles[l].v.x();
				new_vy = particles[l].v.y();
				ghost = Particle(new_x, new_y, new_vx, new_vy, particles[l].rho, particles[l].p, particles[l].m, Ghost);
				addParticleToGrid(ghost, particles.size());
				particles.push_back(ghost);
			}
			for (k = 0; k < grid.elem[grid.wid + res - 1 - i][j].size(); k++)
			{
				l = grid.elem[grid.wid + res - 1 - i][j][k];
				new_x = 2.0 * dom_size - particles[l].x.x();
				// cout << new_x << endl;
				new_y = particles[l].x.y();
				new_vx = - particles[l].v.x();
				new_vy = particles[l].v.y();
				ghost = Particle(new_x, new_y, new_vx, new_vy, particles[l].rho, particles[l].p, particles[l].m, Ghost);
				addParticleToGrid(ghost, particles.size());
				particles.push_back(ghost);
			}
		}
	/*
	for (j = grid.num + 2 * grid.wid - 1; j >= 0; j--)
	{
		for (i = 0; i < grid.num + 2 * grid.wid; i++)
			if (grid.elem[i][j].size() > 0)
				for (int k = 0; k < grid.elem[i][j].size(); k++)
					cout << grid.elem[i][j][k] << ' ';
			else
				cout << ' ';
		cout << endl;
	}
	*/
	
}

void ParticleSystem::deleteGhostParticles()
{
	int i, j;
	for (i = 0; i < 2 * grid.wid + res; i++)
		for (j = 0; j < grid.wid; j++)
		{
			grid.elem[i][j].clear();
			grid.elem[i][grid.wid + res + j].clear();
		}
	for (i = 0; i < grid.wid; i++)
		for (j = grid.wid; j < grid.wid + res; j++)
		{
			grid.elem[i][j].clear();
			grid.elem[grid.wid + res + i][j].clear();
		}
	while (particles.size() > n_real)
		particles.pop_back();
}

void ParticleSystem::handleCollisions()
{
	for (int i = 0; i < n_real; i++)
	{
		if (particles[i].x.x() < 0)
			particles[i].x = Vector2d(0.001, particles[i].x.y());
		if (particles[i].x.x() > 1)
			particles[i].x = Vector2d(0.999, particles[i].x.y());
		if (particles[i].x.y() < 0)
			particles[i].x = Vector2d(particles[i].x.x(), 0.001);
		if (particles[i].x.y() > 1)
			particles[i].x = Vector2d(particles[i].x.x(), 0.999);
	}
}