#include "Grid.h"

Grid::Grid() 
{
}

Grid::~Grid()
{
}

Grid::Grid(int n, int d)
{
	num = n;
	wid = d;
	// allocate space for elem
	elem = new vector<int>* [n + 2 * d];
	for (int i = 0; i < n + 2 * d; i++)
	{
		elem[i] = new vector<int>[n + 2 * d];
		for (int j = 0; j < n + 2 * d; j++)
			elem[i][j].clear();
	}
};