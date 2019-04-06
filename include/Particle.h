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

	Matrix2d		Cp;					// gradient of velocity
	Matrix2d		dg;					// deformation gradient
	double			jc;					// Jacobian

public:
	Particle();
	~Particle();
	Particle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume);
};