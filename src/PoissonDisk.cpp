#include <iostream>
#include <ctime>
#include "PoissonDisk.h"
#include "Parameters.h"

PoissonDisk::PoissonDisk() 
{
}

PoissonDisk::~PoissonDisk()
{
}

PointList PoissonDisk::sample(int trial_num, double min_distance)
{
	int i=0, j=0, p=0, k=0, index=0, new_i=0, new_j=0, new_index=0, ref_i=0, ref_j=0;
	bool flag_1=false, flag_2=false;
	double r=0, theta=0;
	Vector2d x0(0,0), prev_x(0,0), new_x(0,0), dis(0,0);

	srand(time(NULL));
	/* Step 0: initialize grid and active list */
	int cell_num = ceil(sqrt(2) / min_distance);
	double width = 1 / double(cell_num);
	grid = (int*)malloc(cell_num * cell_num * sizeof(int));
	for (i = 0; i < cell_num * cell_num; i++)
		grid[i] = -1;
	PointList sample_particles;
	sample_particles.clear();
	active_list.clear();

	/* Step 1: select initial sample */
	x0 << 0.5, 0.5;
	sample_particles.push_back(x0);
	i = floor(x0(0) / width);
	j = floor(x0(1) / width);
	grid[i * cell_num + j] = 0;
	active_list.push_back(0);

	/* Step 2: sample more particles */
	while (!active_list.empty()) {
		// randomly choose an index
		p = rand() % active_list.size();
		index = active_list[p];
		prev_x = sample_particles[index];
		flag_1 = false;
		for (k = 1; k <= trial_num; k++) {
			// choose particle in radical annulus
			r = min_distance + min_distance * rand() / (double)RAND_MAX;
			theta = 2 * Pi * rand() / (double)RAND_MAX;
			new_x << prev_x(0) + r * cos(theta), prev_x(1) + r * sin(theta);
			// identify grid node
			new_i = floor(new_x(0) / width);
			new_j = floor(new_x(1) / width);
			flag_2 = true;
			if (new_i < 0 || new_i >= cell_num || new_j < 0 || new_j >= cell_num)
				continue;
			if (grid[new_i * cell_num + new_j] >= 0)
				continue;
			for (ref_i = new_i - 1; ref_i <= new_i + 1; ref_i++)
				for (ref_j = new_j - 1; ref_j <= new_j + 1; ref_j++)
				{
					index = grid[ref_i * cell_num + ref_j];
					if (ref_i >= 0 && ref_i < cell_num && ref_j >= 0 && ref_j < cell_num && index >= 0)
					{
						Vector2d ref_pt = sample_particles[index];
						dis = new_x - ref_pt;
						if (dis.norm() < min_distance)
							flag_2 = false;
					}
				}
			if (!flag_2)
				continue;
			// rejection for neighbor particles
			// sample
			sample_particles.push_back(new_x);
			new_index = sample_particles.size() - 1;
			grid[new_i * cell_num + new_j] = new_index;
			active_list.push_back(new_index);
			flag_1 = true;
		}
		// remove p from active list
		if (!flag_1)
			active_list.erase(active_list.begin() + p);
	}
	free(grid);
	return sample_particles;
}