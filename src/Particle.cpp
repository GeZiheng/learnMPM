#include "Particle.h"

using namespace std;
using namespace Eigen;

Particle::Particle() 
{
}

Particle::~Particle() 
{
}

Particle::Particle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	V = volume;
	Cp = MatrixXd::Zero(2, 2);
	dg = MatrixXd::Identity(2, 2);
	jc = 1;
}