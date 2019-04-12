#include "Functions.h"
#include <iostream>
#include <cmath>

using namespace std;

void calWeight(const Vector2d &x, double h, int base_node[2], double weight[3][3])		// calculate weight on 3 by 3 stencil
{
	double wx[3], wy[3];
	base_node[0] = floor(x(0) / h - 0.5);
	base_node[1] = floor(x(1) / h - 0.5);
	double rx = x(0) / h - base_node[0];
	double ry = x(1) / h - base_node[1];
	wx[0] = 0.5 * (1.5 - rx) * (1.5 - rx);
	wx[1] = 0.75 - (rx - 1) * (rx - 1);
	wx[2] = 0.5 * (rx - 0.5) * (rx - 0.5);
	wy[0] = 0.5 * (1.5 - ry) * (1.5 - ry);
	wy[1] = 0.75 - (ry - 1) * (ry - 1);
	wy[2] = 0.5 * (ry - 0.5) * (ry - 0.5);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			weight[i][j] = wx[i] * wy[j];
}

void calWeightGrad(const Vector2d &x, double h, int base_node[2], Vector2d weight_grad[3][3])		// calculate weight gradient on 3 by 3 stencil
{
	double wx[3], dwx[3], wy[3], dwy[3];
	base_node[0] = floor(x(0) / h - 0.5);
	base_node[1] = floor(x(1) / h - 0.5);
	double rx = x(0) / h - base_node[0];
	double ry = x(1) / h - base_node[1];
	wx[0] = 0.5 * (1.5 - rx) * (1.5 - rx);
	wx[1] = 0.75 - (rx - 1) * (rx - 1);
	wx[2] = 0.5 * (rx - 0.5) * (rx - 0.5);
	dwx[0] = (rx - 1.5) / h;
	dwx[1] = 2 * (1 - rx) / h;
	dwx[2] = (rx - 0.5) / h;
	wy[0] = 0.5 * (1.5 - ry) * (1.5 - ry);
	wy[1] = 0.75 - (ry - 1) * (ry - 1);
	wy[2] = 0.5 * (ry - 0.5) * (ry - 0.5);
	dwy[0] = (ry - 1.5) / h;
	dwy[1] = 2 * (1 - ry) / h;
	dwy[2] = (ry - 0.5) / h;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			weight_grad[i][j] << dwx[i] * wy[j], wx[i] * dwy[j];
}

void polarDecompos(const Matrix2d &F, Matrix2d &R, Matrix2d &S)
{
	double e = F(1, 0) - F(0, 1);
	double f = F(0, 0) + F(1, 1);
	double scale = 1 / sqrt(e * e + f * f);
	double c = f * scale;
	double s = e * scale;
	R(0, 0) = c;
	R(0, 1) = -s;
	R(1, 0) = s;
	R(1, 1) = c;

	S(0, 0) = c * F(0, 0) + s * F(1, 0);
	S(0, 1) = c * F(0, 1) + s * F(1, 1);
	S(1, 0) = S(0, 1);
	S(1, 1) = -s * F(0, 1) + c * F(1, 1);
}