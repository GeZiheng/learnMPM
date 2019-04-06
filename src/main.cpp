#include <iostream>
#include <fstream>
#include "Simulator.h"

using namespace std;

int main(int argc, char **argv)
{
	Simulator sim;
	sim.init();
	while (sim.getTime() < sim.end_t)
		sim.oneTimeStep();
}