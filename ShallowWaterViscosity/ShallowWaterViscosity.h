#ifndef _SHALLOW_WATER_VISCOSITY_H_
#define _SHALLOW_WATER_VISCOSITY_H_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<omp.h>
#include <ctime>

#define MACHEPS 1.11022e-16
#define M_PI 3.14159265358979323846
#define gamma 2.

struct Values
{
	double u;
	double rho; /* rho = h */
	double p; /* p = h*h/2/Fr/Fr */
};

struct Values1
{
	double u;
	double v;
	double h;
};

typedef struct Point{
	double x, y;
	double u, v, h;
	double u1, u2, u3;
	char c, square;
} Point;

typedef struct All_data {
	unsigned int n; /*horisontal cells*/
	unsigned int m; /*vertical cells*/
	double d_x, d_y; /*cell edges*/
	double S; /*cell area*/
	double cur_time, d_t;
	double Q_l, Q_r;
	double Fr; /*Frud number*/
	Point** layer;
} All_data;

#endif /* _SHALLOW_WATER_VISCOSITY_H_ */