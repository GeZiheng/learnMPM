#pragma once

#include <vector>
#include <Eigen>

using namespace std;
using namespace Eigen;

class Particle
{
public:
	Vector2d		x;					// position
	Vector2d		v;					// velocity
	double			m;					// mass
	double			V;					// volume
	Matrix2d		v_grad;				// gradient of velocity
	double			v_div;				// divergence of velocity
	Matrix2d		dg;					// deformation gradient
	double			jc;					// Jacobian

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
	Particle();
	~Particle();
	virtual void evolveDG(double dt);
	virtual void calStress(double mu, double lambda, Matrix2d &tmp);
};

class JellyParticle : public Particle
{
public:
	JellyParticle();
	~JellyParticle();
	JellyParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume);
	void evolveDG(double dt);
	void calStress(double mu, double lambda, Matrix2d &tmp);
};

class WaterParticle : public Particle
{
public:
	WaterParticle();
	~WaterParticle();
	WaterParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume);
	void evolveDG(double dt);
	void calStress(double mu, double lambda, Matrix2d &tmp);
};