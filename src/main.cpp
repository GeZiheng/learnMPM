#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include "Parameters.h"
#include "Simulator.h"

using namespace std;

Simulator sim;

void InitGL()
{
	glClearColor(1.0, 1.0, 1.0, 0.9);
	glColor3f(0.0f, 0.0f, 1.0f);
	glPointSize(4.0);
	glEnable(GL_POINT_SMOOTH);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, window_width, 0.0, window_height);
}

void InitSPH()
{
	sim.init();
}

void Update()
{
	sim.oneTimeStep();
	glutPostRedisplay();
	if (sim.getTime() >= end_t)
		exit(0);
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(1.2 * window_width, 1.2 * window_height);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("SPH simulation");
	InitGL();
	InitSPH();
	glutDisplayFunc(Update);
	glutMainLoop();
	return 0;
}