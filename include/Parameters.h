#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <cmath>

const double Pi = 4.0 * atan(1);

extern int res;						// grid res
extern double dx;					// grid dx width
extern double h;					// particle search radius
extern double dt;					// time step length
extern double rho0;					// initial rho0
extern double gamma;				// gamma
extern double gravity;				// gravity
extern double c;					// speed of sound
extern double supp_size;			// support size
extern double water_width;			// width of water
extern double water_height;			// height of water
extern double end_t;				// end time
extern int frame_num;				// number of frame
extern double frame_dt;				// frame dt
extern double window_width;			// width of window
extern double window_height;		// height of window

#endif