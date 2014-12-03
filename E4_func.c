/*
E4_func.c
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

double acceleration(double v, double x, double ny, double omega){

	double a;
	double eta;

	eta = (double) rand() / RAND_MAX;

	a = -ny * v - omega*omega*x + eta;

	return a;

}





