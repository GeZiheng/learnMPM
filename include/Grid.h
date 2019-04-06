#pragma once
#include <iostream>
#include <vector>
#include <Eigen>

using namespace std;
using namespace Eigen;

class GridNode
{
public:
	Vector2d		x;					// position
	Vector2d		v;					// velocity
	double			m;					// mass
	Vector2d		f;					// force
	Vector2d		mv;					// momentum

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
	GridNode();
	~GridNode();
	void SetVal(double pos_x, double pos_y, double velo_x, double velo_y, double mass);
};

class Grid
{
public:
	int			res;			// grid resolution

private:
	double		dx;				// grid dx
	GridNode	*elem;			// storage of grid nodes

public:
	Grid();
	~Grid();
	Grid(int n);
	GridNode* GetElem(int i, int j);		// get grid elem (i,j) pointer
	void Reset();					// reset grid data
};