#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Parameters.h"
#include "Simulator.h"
#include "Functions.h"
#include "PoissonDisk.h"

Simulator::Simulator()
{
}

Simulator::~Simulator()
{
	delete grid;
}

void Simulator::init()
{
	ifstream myfile;
	string s, str;
	int i;
	myfile.open("../config.txt");

	// get rho0 from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	rho0 = stod(str);

	// get young's modulus and poisson's ratio from file
	str.clear();
	getline(myfile, s);
	for (i = 2; i < s.length(); i++)
		str.push_back(s[i]);
	double E = stod(str);
	str.clear();
	getline(myfile, s);
	for (i = 3; i < s.length(); i++)
		str.push_back(s[i]);
	double nu = stod(str);
	mu = E / (2 * (1 + nu));
	lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

	// get gravity from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	gravity = stod(str);

	// get center position from file
	double rx, ry;
	str.clear();
	getline(myfile, s);
	i = 7;
	while (s[i] != ' ')
	{
		str.push_back(s[i]);
		i++;
	}
	rx = stod(str);
	i++;
	str.clear();
	while (i < s.length())
	{
		str.push_back(s[i]);
		i++;
	}
	ry = stod(str);
	rec_center = Vector2d(rx, ry);

	// get rect size from file
	str.clear();
	getline(myfile, s);
	for (i = 5; i < s.length(); i++)
		str.push_back(s[i]);
	rec_size = stod(str);

	// get dx from file
	str.clear();
	getline(myfile, s);
	for (i = 11; i < s.length(); i++)
		str.push_back(s[i]);
	int res = stoi(str);
	dx = 1.0 / res;

	// get frame dt from file
	str.clear();
	getline(myfile, s);
	for (i = 11; i < s.length(); i++)
		str.push_back(s[i]);
	int frame_rate = stoi(str);
	frame_num = 0;
	frame_dt = 1.0 / double(frame_rate);

	// get dt from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	dt = stod(str);

	// get end_t from file
	str.clear();
	getline(myfile, s);
	for (i = 6; i < s.length(); i++)
		str.push_back(s[i]);
	end_t = stod(str);

	myfile.close();

	t = 0;

	// init grid
	grid = new Grid(res);
	// init particles
	PoissonDisk pd;
	vector<Vector2d> pts_pos = pd.sample(30, 0.25 / (res*rec_size));
	pts_num = pts_pos.size();
	for (i = 0; i < pts_num; i++)
	{
		rx = rec_center.x() + rec_size * (pts_pos[i].x() - 0.5);
		ry = rec_center.y() + rec_size * (pts_pos[i].y() - 0.5);
		pts_cloud.push_back(Particle(rx, ry, 0, 0, rho0 * dx * dx / 8, dx * dx / 8));
	}
}

void Simulator::writeData(int frame_num)
{
	ofstream myfile;
	ostringstream filename;
	Vector2d x_i;
	filename << "../output/frame" << setfill('0') << setw(4) << frame_num << ".dat";
	myfile.open(filename.str());
	myfile << "variables=\"x\",\"y\",\"u\",\"v\",\"m\",\"jc\"" << endl;
	for (int i = 0; i < pts_num; i++)
		myfile << scientific << setprecision(16) << pts_cloud[i].x.x() << ' ' << pts_cloud[i].x.y() << ' ' << pts_cloud[i].v.x() << ' ' << pts_cloud[i].v.y() << ' ' << pts_cloud[i].m << ' ' << pts_cloud[i].jc << endl;
	myfile.close();
}

void Simulator::P2G()
{
	int n, bn_i, bn_j, i, j;
	double weight;
	// reset grid values
	grid->Reset();
	for (n = 0; n < pts_num; n++)
	{
		bn_i = floor(pts_cloud[n].x.x() / dx - 0.5);
		bn_j = floor(pts_cloud[n].x.y() / dx - 0.5);
		for (i = bn_i; i < bn_i + 3; i++)
			for(j = bn_j; j < bn_j + 3; j++)
				if (i >= 0 && i <= grid->res && j >= 0 && j <= grid->res)
				{
					weight = calWeight(pts_cloud[n].x - grid->GetElem(i,j)->x, dx);
					grid->GetElem(i,j)->m += pts_cloud[n].m * weight;
					Vector2d pv = pts_cloud[n].v + pts_cloud[n].Cp * (grid->GetElem(i, j)->x - pts_cloud[n].x);
					grid->GetElem(i,j)->mv += pts_cloud[n].m * pv * weight;
				}
	}
	for (i = 0; i <= grid->res; i++)
		for (j = 0; j <= grid->res; j++)
			if (grid->GetElem(i,j)->m > 0)
				grid->GetElem(i,j)->v = grid->GetElem(i,j)->mv / grid->GetElem(i,j)->m;
}

void Simulator::G2P()
{
	int n, bn_i, bn_j, i, j;
	double weight;
	for (n = 0; n < pts_num; n++)
	{
		pts_cloud[n].v.setZero();
		pts_cloud[n].Cp.setZero();
		bn_i = floor(pts_cloud[n].x.x() / dx - 0.5);
		bn_j = floor(pts_cloud[n].x.y() / dx - 0.5);
		for (i = bn_i; i < bn_i + 3; i++)
			for (j = bn_j; j < bn_j + 3; j++)
				if (i >= 0 && i <= grid->res && j >= 0 && j <= grid->res)
				{
					Vector2d dis = pts_cloud[n].x - grid->GetElem(i, j)->x;
					weight = calWeight(pts_cloud[n].x - grid->GetElem(i,j)->x, dx);
					pts_cloud[n].v += grid->GetElem(i,j)->v * weight;
					pts_cloud[n].Cp += grid->GetElem(i,j)->v * (grid->GetElem(i,j)->x - pts_cloud[n].x).transpose() * weight;
				}
		pts_cloud[n].Cp *= 4 / (dx * dx);
	}
}

void Simulator::gridF()
{
	Matrix2d F, U, V, R, P, tmp;
	double J;
	int n, bn_i, bn_j, i, j;
	Vector3d temp;
	Vector2d weight_grad;
	for (n = 0; n < pts_num; n++)
	{
		F = pts_cloud[n].dg;
		J = pts_cloud[n].jc;
		JacobiSVD<Matrix2d> svd(F, ComputeFullU | ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		R = U * V.transpose();
		P = 2 * mu * (F - R) + lambda * J * (J - 1) * F.transpose().inverse();
		tmp = pts_cloud[n].V * P * F.transpose();
		bn_i = floor(pts_cloud[n].x.x() / dx - 0.5);
		bn_j = floor(pts_cloud[n].x.y() / dx - 0.5);
		for (i = bn_i; i < bn_i + 3; i++)
			for (j = bn_j; j < bn_j + 3; j++)
				if (i >= 0 && i <= grid->res && j >= 0 && j <= grid->res)
				{
					temp = calWeightGrad(pts_cloud[n].x - grid->GetElem(i, j)->x, dx);
					weight_grad = Vector2d(temp(1), temp(2));
					grid->GetElem(i, j)->f -= tmp * weight_grad;
				}
	}
}

void Simulator::gridV()
{
	double rx, ry, vx, vy;
	for (int i = 0; i <= grid->res; i++)
		for (int j = 0; j <= grid->res; j++)
			if (grid->GetElem(i, j)->m > 0)
			{
				grid->GetElem(i, j)->v += (grid->GetElem(i,j)->f / grid->GetElem(i,j)->m + Vector2d(0, -gravity)) * dt;
				rx = grid->GetElem(i, j)->x.x();
				ry = grid->GetElem(i, j)->x.y();
				vx = grid->GetElem(i, j)->v.x();
				vy = grid->GetElem(i, j)->v.y();
				// grid collision handling
				if (rx < 2 * dx && vx < 0 || rx > 1 - 2 * dx && vx > 0)
					grid->GetElem(i, j)->v << 0.0, vy;
				if (ry < 2 * dx && vy < 0 || ry > 1 - 2 * dx && vy > 0)
					grid->GetElem(i, j)->v << vx, 0.0;
			}
}

void Simulator::ptsDG()
{
	for (int n = 0; n < pts_num; n++)
	{
		pts_cloud[n].dg *= (Matrix2d::Identity() + pts_cloud[n].Cp * dt);
		pts_cloud[n].jc = pts_cloud[n].dg.determinant();
	}
}

void Simulator::advection()
{
	for (int n = 0; n < pts_num; n++)
		pts_cloud[n].x += pts_cloud[n].v * dt;
}

void Simulator::oneTimeStep()
{
	double tmp_dt = dt;
	if (t + dt > frame_num * frame_dt)
	{
		writeData(frame_num);
		if (t < frame_num * frame_dt)
			dt = frame_num * frame_dt - t;
		frame_num++;
		cout << t << endl;
	}
	P2G();
	gridF();
	gridV();
	ptsDG();
	G2P();
	advection();
	t += dt;
	dt = tmp_dt;
}

double Simulator::getTime()
{
	return t;
}