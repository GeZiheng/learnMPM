#include "Grid.h"

GridNode::GridNode()
{
}

GridNode::~GridNode()
{
}

void GridNode::SetVal(double pos_x, double pos_y, double velo_x, double velo_y, double mass)
{
	x = Vector2d(pos_x, pos_y);
	v = Vector2d(velo_x, velo_y);
	m = mass;
	f = Vector2d(0, 0);
	mv = m * v;
}


Grid::Grid() 
{
}

Grid::~Grid()
{
	delete[] elem;
}

Grid::Grid(int n)
{
	res = n;
	dx = 1.0 / n;
	// allocate space for elem
	elem = new GridNode [(n+1)*(n+1)];
	int i, j;
	for (i = 0; i <= n; i++)
		for (j = 0; j <= n; j++)
			GetElem(i, j)->SetVal(1.0 * i / n, 1.0 * j / n, 0.0, 0.0, 0.0);
}

void Grid::Reset()
{
	for (int i = 0; i <= res; i++)
		for (int j = 0; j <= res; j++)
		{
			GetElem(i, j)->m = 0.0;
			GetElem(i, j)->v.setZero();
			GetElem(i, j)->mv.setZero();
			GetElem(i, j)->f.setZero();
		}
}

GridNode* Grid::GetElem(int i, int j)
{
	return &elem[i*(res+1)+j];
}