#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Grid
{
public:
	int			num;			// number of grid nodes in domain
	int			wid;			// number of grid nodes in each side of extended domain

public:
	vector<int> **elem;

public:
	Grid();
	~Grid();
	Grid(int n, int d);
};