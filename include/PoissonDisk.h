#pragma once

#include <vector>
#include <Eigen>

using namespace std;
using namespace Eigen;

typedef vector<Vector2d, aligned_allocator<Vector2d>> PointList;

class PoissonDisk 
{
public:
	int* grid;
	vector<int> active_list;

public:
	PoissonDisk();
	~PoissonDisk();
	PointList sample(int trial_num, double min_distance);
};