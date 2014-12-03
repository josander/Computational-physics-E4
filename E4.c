/*
 E4.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "E4_func.h"
#define PI 3.141592653589
# define kB 0.000086173324

/* Main program */
int main()
{

	// Declaration of variables
	int i, j, k;
	double ny;
	double omega;
	double c0;
	double m;
	double temp;
	double dt;
	double timesteps;
	double g1, g2;

	// Initiation of variables
	ny = 1;
	omega = 1;
	c0 = exp(-ny * dt);
	m = 1;
	temp = 300;
	dt = 0.1;
	timesteps = 1000;
	
	// Declaration of arrays
	double *x = malloc(timesteps * sizeof(double));
	double *v = malloc(timesteps * sizeof(double));
	double *a = malloc(timesteps * sizeof(double));

	// Initiation of arrays
	for(j = 0; j < timesteps; j++){
		x[j] = 0.0;
		v[j] = 0.0;
		a[j] = 0.0;
	}

	// Seed for generating random numbers
	srand(time(NULL));

	// Time evolution
	for(i = 0; i < timesteps; i++){

		// Generate two random numbers with gaussian distribution
		g1 = (double) rand() / RAND_MAX;
		g2 = (double) rand() / RAND_MAX;		
	
		// v(t+)
		v[i] = sqrt(c0)*v[i] + sqrt(kB * temp / m) * sqrt(1-c0)*g1;
	
		// v(t + dt/2)
		v[i] = v[i] + a[i] * dt / 2;
	
		// x(t + dt)
		x[i] = x[i] + v[i] * dt;

		// Calc new acceleration
		a[i] = acceleration(v[i], x[i], ny, omega);

		// v(t- + dt)
		v[i] = v[i] + a[i] * dt / 2;

		// v(t + dt)
		v[i] = sqrt(c0) * v[i] + sqrt(kB * temp / m) * sqrt(1-c0)*g2;

	}

	// Free allocated memory
	free(x); free(v); free(a);
	x = NULL; v = NULL; a = NULL;

}

