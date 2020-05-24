#include "Particle.h"
#include "Functions.h"
#include <iostream>

using namespace std;
using namespace Eigen;

Particle::Particle() 
{
}

Particle::~Particle() 
{
}

void Particle::evolveDG(double dt)
{
}

void Particle::calStress(double mu, double lambda, Matrix2d &tmp)
{
}

JellyParticle::JellyParticle()
{
}

JellyParticle::JellyParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	V = volume;
	v_grad = MatrixXd::Zero(2, 2);
	v_div = 0;
	dg = MatrixXd::Identity(2, 2);
	jc = 1;
}


JellyParticle::~JellyParticle()
{
}

void JellyParticle::calStress(double mu, double lambda, Matrix2d & tmp)
{
	Matrix2d R, S, P;
	R.setZero();
	polarDecompos(dg, R, S);
	// JacobiSVD<Matrix2d> svd(dg, ComputeFullU | ComputeFullV);
	// R = svd.matrixU() * svd.matrixV().transpose();
	P = 2 * mu * (dg - R) + lambda * jc * (jc - 1) * dg.transpose().inverse();
	tmp = V * P * dg.transpose();
}

void JellyParticle::evolveDG(double dt)
{
	dg = (Matrix2d::Identity() + v_grad * dt) * dg;
	jc *= 1 + v_div * dt;
}

WaterParticle::WaterParticle()
{
}

WaterParticle::WaterParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	V = volume;
	v_grad = MatrixXd::Zero(2, 2);
	v_div = 0;
	dg = MatrixXd::Identity(2, 2);
	jc = 1;
}


WaterParticle::~WaterParticle()
{
}

void WaterParticle::calStress(double mu, double lambda, Matrix2d &tmp)
{
	Matrix2d sigma = mu * (1 - pow(jc, -lambda)) * Matrix2d::Identity();
	tmp = V * jc * sigma;
}

void WaterParticle::evolveDG(double dt)
{
	jc *= 1 + v_div * dt;
}

SnowParticle::SnowParticle()
{
}

SnowParticle::SnowParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	V = volume;
	v_grad = MatrixXd::Zero(2, 2);
	v_div = 0;
	dg = MatrixXd::Identity(2, 2);
	jc = 1;
}


SnowParticle::~SnowParticle()
{
}

void SnowParticle::calStress(double mu, double lambda, Matrix2d &tmp)
{
	Matrix2d R, S, P;
	R.setZero();
	JacobiSVD<Matrix2d> svd(dg, ComputeFullU | ComputeFullV);
	R = svd.matrixU() * svd.matrixV().transpose();
	P = 2 * mu * (dg - R) + lambda * jc * (jc - 1) * dg.transpose().inverse();
	tmp = V * P * dg.transpose();
}

void SnowParticle::evolveDG(double dt)
{
	dg = (Matrix2d::Identity() + v_grad * dt) * dg;
	jc *= 1 + v_div * dt;
}

SandParticle::SandParticle()
{
}

SandParticle::SandParticle(double pos_x, double pos_y, double velo_x, double velo_y, double mass, double volume)
>>>>>>> Fixed code
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	V = volume;
	v_grad = MatrixXd::Zero(2, 2);
	v_div = 0;
	dg = MatrixXd::Identity(2, 2);
	jc = 1;
}

SandParticle::~SandParticle()
{
}

void SandParticle::calStress(double mu, double lambda, Matrix2d &tmp)
{
	Matrix2d R, S, P;
	R.setZero();
	JacobiSVD<Matrix2d> svd(dg, ComputeFullU | ComputeFullV);
	R = svd.matrixU() * svd.matrixV().transpose();
	P = 2 * mu * (dg - R) + lambda * jc * (jc - 1) * dg.transpose().inverse();
	tmp = V * P * dg.transpose();
}

void SandParticle::evolveDG(double dt)
{
	dg = (Matrix2d::Identity() + v_grad * dt) * dg;
	jc *= 1 + v_div * dt;
}