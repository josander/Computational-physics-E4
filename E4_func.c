/*
E4_func.c
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

// Function that evaluates the integral of x(x+1) with the Monte Carlo method. This method uses a variable with a uniform distribution between 0 and 1.
void integral_uniform(){
	
	// Declaration and initiation of variables
	int i, j;
	double sum, sum2;
	int N = 1;
	double mean, mean2;
	double var;

	// Print the expected value of the integral
	printf("UNIFORM: Expected value: %.8F ± %.8F \n", (double) 1/6, 0.03726852);

	// For N = 10, 100, 1000, 10000
	for(i = 0; i < 4; i++){

		N = N*10;
		double x[N];
		sum = 0;
		sum2 = 0;
		var = 0;
		
		// Calculate the integral
		for(j = 0; j < N; j++){

			x[j] = ((double) rand()/ (double) RAND_MAX);
			sum += x[j]*(1-x[j]);
			sum2 += x[j]*(1-x[j]) * x[j]*(1-x[j]);

		}

	// Get the mean
	mean = sum/N;
	mean2 = sum2/N;
	var = (mean2 - mean*mean)/N;	

	// Print the result in the terminal
	printf("For N = %i \t Integral = %.8F ± %.8F \n", N, mean, sqrt(var));

	}

}

// Function that evaluates the integral of x(x+1) with the Monte Carlo method. This function uses importance sampling where the variable has a sine-distribution between 0 and 1.
void integral_sine(){
	
	// Declaration and initiation of variables
	int i, j;
	double sum, sum2;
	int N = 1;
	double mean, mean2;
	double var;
	double x[4][10000];

	// Open a file to print the variable x in
	FILE *x_file;
	x_file = fopen("distribution.data","w");

	// Print the expected value of the integral
	printf("SINE: Expected: %.8F ± %.8F \n", (double) 1/6, 0.03726852);

	// For N = 10, 100, 1000, 10000
	for(i = 0; i < 4; i++){

		N = N*10;
		sum = 0;
		sum2 = 0;
		var = 0;
		
		// Calculate the integral
		for(j = 0; j < N; j++){

			// Random numbers with  a sinusiodal distribution
			x[i][j] = ((double) rand()/ (double) RAND_MAX);
			x[i][j] = acos(1 - 2 * x[i][j])/PI;

			sum += x[i][j] * (1-x[i][j]) * 2 / sin(PI * x[i][j]) / PI;
			sum2 += x[i][j]*(1-x[i][j]) * 2 / sin(PI*x[i][j]) * x[i][j]*(1-x[i][j]) * 2 / sin(PI*x[i][j]) / PI / PI;

		}

		// Get the mean
		mean = sum/N;
		mean2 = sum2/N;
		var = (mean2 - mean*mean)/N;	

		// Print the result in the terminal
		printf("For N = %i \t Integral = %.8F ± %.8F \n", N, mean, sqrt(var));

	}

	// Print x to distribution.data
	for(j = 0; j < N; j++){
		fprintf(x_file,"%F \t %F \t %F \t %F \n", x[0][j], x[1][j], x[2][j], x[3][j]);
	}

	// Close the data-file
	fclose(x_file); 

}

// Function that evaluates the integral of x(x+1) with the Monte Carlo method. This function uses importance sampling where the variable has a sine-distribution between 0 and 1.
void integral_metropolis(){
	
	// Declaration and initiation of variables
	int i, j;
	double sum, sum2;
	int N = 10000;
	double mean, mean2;
	double var;
	double x[N];
	double p[N];
	double delta;
	double q;
	int throw_away, rejections;
	double r;

	// Open a file to print the variable x in
	FILE *m_file;
	m_file = fopen("metropolis.data","w");

	// Print the expected value of the integral
	printf("METROPOLIS: Expected value: %.8F ± %.8F \n", (double) 1/6, 0.03726852);

	// Initialize variables
	sum = 0;
	sum2 = 0;
	var = 0;
	x[0] = 0.5;
	p[0] = sin(PI * x[0]);
	delta = 0.45;
	throw_away = 2000;
	rejections = 0;

	fprintf(m_file,"%F \n", x[0]);
		
	int n = 0;

	// Calculate the integral
	for(j = 1; j < N; j++){

		// Generate random number and get next state x
		r = (double) rand() / (double) RAND_MAX;	
		x[j] = x[j-1] + delta*(r - 0.5);

		// Calculate the probability
		p[j] = 0.0;
		p[j] = sin(PI * x[j]);
		q = p[j]/p[j-1];

		// New random number
		r = (double) rand() / (double) RAND_MAX;

		// Trial
		if(q < r){
			x[j] = x[j-1];
			p[j] = p[j-1];
			rejections++;
		}
		
		// Skip the 'throw_away' first datapoints
		if(j > throw_away){
			n++;
			sum += x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI;
			sum2 += x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI * x[j] * (1-x[j]) * 2.0 / sin(PI * x[j]) / PI;
		}

	}

	// Get the means to calculate the variance
	mean = sum/(N-throw_away-1);
	mean2 = sum2/(N-throw_away-1);
	var = (mean2 - mean*mean)/(N-throw_away-1);	
	
	// In the terminal, print how many throw aways
	printf("Nbr of rejections: %i \n", rejections);

	// Print the result in the terminal
	printf("For N = %i \t Integral = %.8F ± %.8F \n", N-throw_away, mean, sqrt(var));

	// Print x to distribution.data
	for(j = 0; j < N; j++){
		fprintf(m_file,"%F \n", x[j]);
	}

	// Close the data-file
	fclose(m_file); 

}

// Function that caculated the auto-correlation function, the statistical inefficiency
void error_corr_func(float *A, int length){

	// Declaration and initiation of variables
	int i, k;
	double mean = 0;
	double mean2 = 0;
	double s = 0;
	double sigmaTot;
	int steps = 200;

	// Declaration of arrays
	double first_term[steps];
	double corr_func[steps];

	// Initiate the array first_term
	for(k = 0; k < steps; k++){
		first_term[k] = 0.0;
	}

	// Calculate all the expected values of A
	for(i = 0; i < length-steps; i++){
		mean += A[i]/(length-steps);
		mean2 += ((A[i]*A[i])/(length-steps)); 
	}

	// Calculate the first term
	for(i = 0; i < (length-steps); i++){
		for(k = 0; k < steps; k++){
			first_term[k] += (A[i]*A[i+k])/(length-steps);
		}
	}

	// Calculate the correlation function
	for(k = 0; k < steps; k++){
		corr_func[k] = ((first_term[k] - (mean*mean))/(mean2 - (mean*mean)));
	}

	// Calculate the statistical inefficiency
	i = 0;
	while(corr_func[i] >= exp(-2)){
		i++;
	}
	
	s = i;

	sigmaTot = sqrt((mean2 - mean*mean)/steps*s);
	printf("Statistical inefficiency: %F \n", s);
	printf("Result: %.8e ± %.10e \n", mean, sigmaTot);

}


// Calculate the statistical inefficiency from the block average
error_block_average(float *A, int length){

	// Declaration and initiation of variables
	int i, j;
	int block_size;
	double mean, mean2, var_f;
	double mean_F, mean2_F, var_F;
	double s;

	FILE *block;
	block = fopen("block_s.txt","w");

	fprintf(block,"%f \n", 0);
	
	// Calculate statistical inefficiency for different block sizes
	for(block_size = 10; block_size < 1000; block_size = block_size + 10){

		int nbr_blocks = length/block_size;
		double block_means[nbr_blocks];
		double block_means2[nbr_blocks];

		// Determine variance for the whole array
		mean = 0.0;
		mean2 = 0.0;	
		for(i = 0; i < length; i++){
			mean += A[i]/length;
			mean2 += A[i]*A[i]/length;
		}

		var_f = (mean2 - mean*mean);
		//printf("var: %.10f \n", var_f);
	
		// Determine average in each block
		for(i = 0; i < nbr_blocks; i++){
			
			block_means[i] = 0.0;
			block_means2[i] = 0.0;
			
			for(j = 0; j < block_size; j++){
				block_means[i] += A[(i*block_size+j)]/block_size;
			}
		}
	

		// Determine average in each block
		mean_F = 0.0;
		mean2_F = 0.0;
		for(i = 0; i < nbr_blocks; i++){
			mean_F += block_means[i]/nbr_blocks;
			mean2_F += block_means[i]*block_means[i]/nbr_blocks;
		}

		var_F = (mean2_F - mean_F * mean_F);
		s = block_size * var_F / var_f;
		
		fprintf(block,"%f \n", s);

	}

	printf("Statistical inefficiency: %f \n", s);

	fclose(block);

}






