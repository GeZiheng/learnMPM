#pragma once
#include <vector>
#include <Eigen>

using namespace std;
using namespace Eigen;

enum ParticleType
{
	Real,
	Ghost
};

class Particle
{
public:
	Vector2d		x;					// position
	Vector2d		v;					// velocity
	double			rho;				// rho0
	double			p;					// pressure
	double			m;					// mass
	ParticleType	type;				// particle type: real or ghost
	int				ref_index;			// reference index (for ghost particle)

public:
	Particle();
	~Particle();
	Particle(double pos_x, double pos_y, double velo_x, double velo_y, double rho1, double p1, double mass, ParticleType t, int i);
};