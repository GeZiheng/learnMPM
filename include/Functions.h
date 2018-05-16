#pragma once
#include "Parameters.h"
#include <iostream>
#include "Particle.h"

double weightFunction(Vector2d x)					// weight function
{
	double s = x.norm() / h;
	double c = 1.0 / ( (Pi * h * h) * (1.0 - 10.0 * exp(-9.0)) );
	if (s < 3)
		return c * (exp(-s * s) - exp(-9.0));
	else
		return 0;
}

Vector2d weightGradient(Vector2d x)				// weight gradient
{
	double r = x.norm();
	double s = r / h;
	double c = 1.0 / ( (Pi * h * h) * (1.0 - 10.0 * exp(-9.0)) );
	if (s > 0 && s < 3)
		return c * (-2 * s * exp(-s * s)) * x / (r*h);
	else
		return Vector2d(0,0);
}

double equation_of_state(double rho)
{
	return rho0 * c * c / gamma * (pow(rho / rho0, gamma) - 1.0);
}

double drhodt_increment(Particle p1, Particle p2)
{
	Vector2d dw = weightGradient(p1.x - p2.x);
	return p2.m * p1.rho / p2.rho * (p1.v - p2.v).dot(dw);
}

double artificial_viscosity(Particle p1, Particle p2)
{
	double phi = 0.1 * h;
	double mu = h * (p1.v - p2.v).dot(p1.x - p2.x) / ((p1.x - p2.x).dot(p1.x - p2.x) + phi * phi);
	double alpha = 0.01;
	double beta = 1.0;
	double rho_ij = (p1.rho + p2.rho) / 2;
	return (-alpha * c * mu + beta * mu * mu) / rho_ij;
}

Vector2d boundary_force(Vector2d x1, Vector2d x2, double r0)
{
	double r = (x1 - x2).norm();
	double q = r0 / r;
	double D = 0.01;
	if (q >= 1)
		return D * (pow(q, 4) - pow(q, 2)) * (x1 - x2) / (r * r);
	else
		return Vector2d(0, 0);
}

Vector2d dvdt_increment(Particle p1, Particle p2)
{
	Vector2d dw = weightGradient(p1.x - p2.x);
	Vector2d incre = -p2.m * (p1.p + p2.p) / (p1.rho * p2.rho) * dw;
	// Vector2d incre = -p2.m * (p1.p / (p1.rho * p1.rho) + p2.p / (p2.rho * p2.rho) ) * dw;
	//if (p2.type == Real)
		//incre += artificial_viscosity(p1, p2) / p1.m;
	if (p2.type == Ghost)
		incre += boundary_force(p1.x, p2.x, dx);
	return incre;
}