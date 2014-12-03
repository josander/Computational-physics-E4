/*
 E4.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "E4_func.h"
#define PI 3.141592653589

/* Main program */
int main()
{
	// Seed for generating random numbers
	srand(time(NULL));

	// Task 1: calculate the integral for different N
	integral_uniform();

	// Task 2: calculate the integral for different N using importance sampling
	integral_sine();

	// Task 3: calculate the integral for different N using the Metropolis algorithm
	integral_metropolis();

	// Task 4: Evaluate the statistical inefficiency
	// Get data from file MC.txt and save it in 'data'
	FILE *fr;
	fr = fopen("MC.txt","r");
	
	int N = 1000000;
    	float data[N];
    	int i;
	
    	for (i = 0; i < N; i++){
		data[i] = 0.0;
        	fscanf(fr, "%f\n", &data[i]);
    	}
	fclose(fr);
	
	// Call the error-method to calculate the statistical inefficiency
	printf("STATISTICAL INEFFICIENCY FROM CORR FUNC\n");	
	error_corr_func(data, N);

	// Call the error-method to calculate the statistical inefficiency
	printf("STATISTICAL INEFFICIENCY FROM BLOCK AVERAGE\n");	
	error_block_average(data, N);

	printf("Done! \n");
}

