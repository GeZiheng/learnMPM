#pragma once
#include <iostream>
#include <Eigen>

using namespace Eigen;

double calWeight(Vector2d x, double h)
{
	double rx, ry, wx, wy;
	rx = x(0) / h;
	ry = x(1) / h;
	if (abs(rx) < 0.5)
		wx = 0.75 - rx * rx;
	else if (abs(rx) < 1.5)
		wx = 0.5 * (1.5 - abs(rx)) * (1.5 - abs(rx));
	else
		wx = 0;
	if (abs(ry) < 0.5)
		wy = 0.75 - ry * ry;
	else if (abs(ry) < 1.5)
		wy = 0.5 * (1.5 - abs(ry)) * (1.5 - abs(ry));
	else
		wy = 0;
	return wx * wy;
}

Vector3d calWeightGrad(Vector2d x, double h)					// weight function
{
	double rx, ry, wx, wy, dwx, dwy;
	rx = x(0) / h;
	if (abs(rx) < 0.5)
	{
		wx = 0.75 - rx * rx;
		dwx = -2 / h * rx;
	}
	else if (abs(rx) < 1.5)
	{
		wx = 0.5 * (1.5 - abs(rx)) * (1.5 - abs(rx));
		dwx = 1 / h * (rx - rx / abs(rx) * 1.5);
	}
	else
	{
		wx = 0;
		dwx = 0;
	}
	ry = x(1) / h;
	if (abs(ry) < 0.5)
	{
		wy = 0.75 - ry * ry;
		dwy = -2 / h * ry;
	}
	else if (abs(ry) < 1.5)
	{
		wy = 0.5 * (1.5 - abs(ry)) * (1.5 - abs(ry));
		dwy = 1 / h * (ry - ry / abs(ry) * 1.5);
	}
	else
	{
		wy = 0;
		dwy = 0;
	}
	return Vector3d(wx * wy, dwx * wy, wx * dwy);
}