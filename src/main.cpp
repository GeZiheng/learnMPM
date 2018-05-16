#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include "Parameters.h"
#include "Simulator.h"

using namespace std;

int main(int argc, char **argv)
{
	Simulator sim;
	sim.init();
	while (sim.getTime() < end_t)
		sim.oneTimeStep();
	return 0;
}