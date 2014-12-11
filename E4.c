/*
 E4.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "E4_func.h"
#define PI 3.141592653589
#define kB (1.3806488 * pow(10,-23)) // Units: [m2 kg /(s2 K)]

/* Main program */
int main()
{

	// Declaration of variables
	int i, j, k;
	double eta;
	double omega;
	double c0;
	double m;
	double temp;
	double dt;
	int timesteps;
	double r1, r2;
	double g1, g2;
	double x, v, a;
	int lags;

	// Initiation of variables
	omega = 3.0;
	eta = 5.0 * omega;
	m = 1.0; // Units: [kg]
	temp = 300.0; // Units: [K]
	dt = 0.05; // Units: [s]
	c0 = exp(-eta * dt);
	timesteps = 100000;
	x = 0.05; // Units: [m]
	v = 0.0;
	a = 0.0;
	double X[timesteps + 1];
	X[0] = x;

	// Allocate memory for vectors
	double *corr_func = malloc((timesteps-500+1) * sizeof(double));

	// Seed for generating random numbers
	srand(time(NULL));

	// File to print distributions
	FILE *dist;
	dist = fopen("distribution.data","w");	

	// File to print the trajectory
	FILE *tr;
	tr = fopen("trajectory5.data","w");	

	// File to print the trajectory
	FILE *corr;
	corr = fopen("corrfunc5.data","w");	

	// Print the trajectory
	fprintf(tr, "%f \n", x);

	// Calc initial acceleration
	a = acceleration(x, omega);

	// Time evolution
	for(i = 1; i < timesteps + 1; i++){

		// Generate two random numbers with uniform distribution
		r1 = (double) rand() / RAND_MAX;
		r2 = (double) rand() / RAND_MAX;

		// Bow-Muller transformation to get random variables with standard normal distribution
		g1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
		g2 = sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);

		// Print all random variables
		fprintf(dist, "%f \t %f \t %f \t %f \n", r1, r2, g1, g2);
	
		// v(t+)
		v = sqrt(c0)*v + sqrt(kB * temp / m) * sqrt(1.0-c0)*g1;
	
		// v(t + dt/2)
		v = v + a * dt / 2.0;
	
		// x(t + dt)
		x = x + v * dt;

		// Calc new acceleration
		a = acceleration(x, omega);

		// v(t- + dt)
		v = v + a * dt / 2.0;

		// v(t + dt)
		v = sqrt(c0) * v + sqrt(kB * temp / m) * sqrt(1-c0)*g2;

		// Print the trajectory
		fprintf(tr, "%f \n", x);

		X[i] = x;

	}

	// Calculate the correlation function
	for(i = 0; i < timesteps-500; i++){
		for(k = 0; k < 500; k++){
			corr_func[k] += (X[i]*X[i+k]);
		}
	}

	// Print correlation function to data-file
	for(k = 0; k < 500; k++){
		fprintf(corr,"%f\n", corr_func[k]/(timesteps-500));
	}


	// Close files
	fclose(dist);
	fclose(tr);
	fclose(corr);

}

