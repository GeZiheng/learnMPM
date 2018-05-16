#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Parameters.h"
#include "Simulator.h"
#include "Functions.h"

// global variables
int res;					// grid res
double dx;					// grid dx width
double h;					// particle search radius
double dt;					// time step length
double rho0;				// initial rho0
double gamma;				// gamma
double gravity;				// gravity
double c;					// speed of sound
double supp_size;			// support size
double water_width;			// width of water
double water_height;		// height of water
double dom_size;			// size of simulation domain
double end_t;				// end time
int frame_num;				// number of frame
double frame_dt;			// frame dt

Simulator::Simulator()
{
}

Simulator::~Simulator()
{
}

void Simulator::init()
{
	ifstream myfile;
	string s, str;
	int i, j;
	myfile.open("../config.txt");

	// get rho0 from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	rho0 = stod(str);

	// get gamma from file
	str.clear();
	getline(myfile, s);
	for (i = 6; i < s.length(); i++)
		str.push_back(s[i]);
	gamma = stod(str);

	// get gravity from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	gravity = stod(str);

	// get support size from file
	str.clear();
	getline(myfile, s);
	for (i = 10; i < s.length(); i++)
		str.push_back(s[i]);
	supp_size = stod(str);

	// get dx from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	dx = stod(str);

	h = supp_size * dx;

	// get domain size from file
	str.clear();
	getline(myfile, s);
	for (i = 9; i < s.length(); i++)
		str.push_back(s[i]);
	dom_size = stod(str);
	res = ceil(dom_size / dx);

	// get water size from file
	str.clear();
	getline(myfile, s);
	for (i = 12; i < s.length(); i++)
		str.push_back(s[i]);
	water_width = stod(str);
	str.clear();
	getline(myfile, s);
	for (i = 13; i < s.length(); i++)
		str.push_back(s[i]);
	water_height = stod(str);

	c = 10 * sqrt(2 * gravity * water_height);

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

	pts_sys = ParticleSystem(res, ceil(3 * supp_size));
	pts_sys_pre = ParticleSystem(res, ceil(3 * supp_size));
	for (j = 0; j < water_height / dx; j++)
		for (i = 0; i < water_width / dx; i++)
		{
			pts_sys.particles.push_back(Particle((i + 0.5) * dx, (j + 0.5) * dx, 0, 0, rho0, 0, rho0*dx*dx, Real));
			pts_sys.n_real ++;
			pts_sys_pre.particles.push_back(Particle((i + 0.5) * dx, (j + 0.5) * dx, 0, 0, rho0, 0, rho0*dx*dx, Real));
			pts_sys_pre.n_real ++;
		}
}

void Simulator::predictor()
{
	vector<int> neighbors;
	double drho_dt;
	Vector2d dv_dt, x_i;
	int i, j;

	for (i = 0; i < pts_sys.n_real; i++)
		pts_sys.particles[i].p = equation_of_state(pts_sys.particles[i].rho);

	pts_sys.countParticlesOnGrid();
	pts_sys.addGhostParticles();

	/* evolve particle rho0 */
	for (i = 0; i < pts_sys.n_real; i++)
	{
		x_i = pts_sys.particles[i].x;
		drho_dt = 0;
		neighbors = pts_sys.findNeighbors(i, 3 * h);
		for (j = 0; j < neighbors.size(); j++)
			drho_dt += drhodt_increment(pts_sys.particles[i], pts_sys.particles[neighbors[j]]);
		pts_sys_pre.particles[i].rho = pts_sys.particles[i].rho + drho_dt * dt / 2;
	}	

	/* evolve particle velocity */
	for (i = 0; i < pts_sys.n_real; i++)
	{
		x_i = pts_sys.particles[i].x;
		dv_dt = Vector2d(0.0, -gravity);
		neighbors = pts_sys.findNeighbors(i, 3 * h);
		for (j = 0; j < neighbors.size(); j++)
			dv_dt += dvdt_increment(pts_sys.particles[i], pts_sys.particles[neighbors[j]]);
		pts_sys_pre.particles[i].v = pts_sys.particles[i].v + dv_dt * dt / 2;
		pts_sys_pre.particles[i].x = pts_sys.particles[i].x + pts_sys.particles[i].v * dt / 2;				// advect particle position
	}

	pts_sys.deleteGhostParticles();
	// pts_sys_pre.handleCollisions();
}

void Simulator::corrector()
{
	vector<int> neighbors;
	double drho_dt;
	Vector2d dv_dt, dw, x_i;
	int i, j;

	for (i = 0; i < pts_sys_pre.n_real; i++)
		pts_sys_pre.particles[i].p = equation_of_state(pts_sys_pre.particles[i].rho);

	pts_sys_pre.countParticlesOnGrid();
	pts_sys_pre.addGhostParticles();

	/* evolve particle rho0 */
	for (i = 0; i < pts_sys.n_real; i++)
	{
		x_i = pts_sys_pre.particles[i].x;
		drho_dt = 0;
		neighbors = pts_sys_pre.findNeighbors(i, 3 * h);
		for (j = 0; j < neighbors.size(); j++)
		{
			// cout << pts_sys_pre.particles[neighbors[j]].x << endl;
			drho_dt += drhodt_increment(pts_sys_pre.particles[i], pts_sys_pre.particles[neighbors[j]]);
		}
		pts_sys.particles[i].rho += drho_dt * dt;
	}

	/* evolve particle velocity */
	for (i = 0; i < pts_sys.n_real; i++)
	{
		x_i = pts_sys_pre.particles[i].x;
		dv_dt = Vector2d(0.0, -gravity);
		neighbors = pts_sys_pre.findNeighbors(i, 3 * h);
		for (j = 0; j < neighbors.size(); j++)
			dv_dt += dvdt_increment(pts_sys_pre.particles[i], pts_sys_pre.particles[neighbors[j]]);
		pts_sys.particles[i].v += dv_dt * dt;
		pts_sys.particles[i].x += pts_sys_pre.particles[i].v * dt;				// advect particle position
	}

	pts_sys_pre.deleteGhostParticles();
	// pts_sys.handleCollisions();

	for (i = 0; i < pts_sys.n_real; i++)
		pts_sys.particles[i].p = equation_of_state(pts_sys.particles[i].rho);
}

void Simulator::writeData(int frame_num)
{
	ofstream myfile;
	ostringstream filename;
	Vector2d x_i;
	filename << "../output/frame" << setfill('0') << setw(4) << frame_num << ".dat";
	myfile.open(filename.str());
	myfile << "variables=\"x\",\"y\",\"p\",\"rho\",\"u\",\"v\"" << endl;
	for (int i = 0; i < pts_sys.n_real; i++)
		myfile << scientific << setprecision(16) << pts_sys.particles[i].x.x() << ' ' << pts_sys.particles[i].x.y() << ' ' << pts_sys.particles[i].p << ' ' << pts_sys.particles[i].rho << ' ' << pts_sys.particles[i].v.x() << ' ' << pts_sys.particles[i].v.y() << endl;
	myfile.close();
}

void Simulator::writeDataa(int frame_num)
{
	ofstream myfile;
	ostringstream filename;
	Vector2d x_i;
	filename << "../output/framee" << setfill('0') << setw(4) << frame_num << ".dat";
	myfile.open(filename.str());
	myfile << "variables=\"x\",\"y\",\"p\",\"rho\",\"u\",\"v\"" << endl;
	for (int i = 0; i < pts_sys.n_real; i++)
		myfile << scientific << setprecision(16) << pts_sys_pre.particles[i].x.x() << ' ' << pts_sys_pre.particles[i].x.y() << ' ' << pts_sys_pre.particles[i].p << ' ' << pts_sys_pre.particles[i].rho << ' ' << pts_sys_pre.particles[i].v.x() << ' ' << pts_sys_pre.particles[i].v.y() << endl;
	myfile.close();
}

void Simulator::densityRecons()
{
	double rho_sum, weight_sum, w;
	int i, j;
	vector<int> neighbors;
	for (i = 0; i < pts_sys.n_real; i++)
	{
		rho_sum = 0;
		weight_sum = 0;
		neighbors = pts_sys.findNeighbors(i, 3.0 * h);
		for (j = 0; j < neighbors.size(); j++)
		{
			w = weightFunction(pts_sys.particles[i].x - pts_sys.particles[neighbors[j]].x);
			rho_sum += pts_sys.particles[neighbors[j]].rho * w;
			weight_sum += w;
		}
		pts_sys.particles[i].rho = rho_sum / weight_sum;
	}

}

void Simulator::oneTimeStep()
{
	double tmp_dt = dt;
	if (t + dt > frame_num * frame_dt)
	{
		writeData(frame_num);
		//writeDataa(frame_num);
		frame_num++;
		//if (t < frame_num * frame_dt)
		//	dt = frame_num * frame_dt - t;
		// densityRecons();
	}
	predictor();
	corrector();
	t += dt;
	dt = tmp_dt;
	cout << t << endl;
}

double Simulator::getTime()
{
	return t;
}