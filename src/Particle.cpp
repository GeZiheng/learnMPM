#include "Particle.h"
#include "Parameters.h"

using namespace std;
using namespace Eigen;

Particle::Particle() 
{
}

Particle::~Particle() 
{
}

Particle::Particle(double pos_x, double pos_y, double velo_x, double velo_y, double rho1, double p1, double mass, ParticleType t)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	rho = rho1;
	p = p1;
	m = mass;
	type = t;
}