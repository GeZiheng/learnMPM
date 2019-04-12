#pragma once

#include <iostream>
#include <Eigen>

using namespace Eigen;

const double Pi = 4.0 * atan(1);

void calWeight(const Vector2d &x, double h, int base_node[2], double weight[3][3]);					// calculate weight on 3 by 3 stencil
void calWeightGrad(const Vector2d &x, double h, int base_node[2], Vector2d weight_grad[3][3]);		// calculate weight gradient on 3 by 3 stencil
void polarDecompos(const Matrix2d &F, Matrix2d &R, Matrix2d &S);											// 2d polar decomposition