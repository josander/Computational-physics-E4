/*	fft.func 	
	Program with (fast) discrete Fourier transform  
	Created by Martin Gren on 2014-10-22
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define PI 3.141592653589

/* Makes fft of data and returns the powerspectrum in fft_data */
void fft(double *data, double *fft_data, int n) /* input data, output fft_data, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	double complex_coefficient[2*n];
	double data_cp[n]; 

	/*make copy of data to avoid messing with data in the transform*/
	for (i = 0; i < n; i++)
	{
		data_cp[i] = data[i];
	}

	/*Declare wavetable and workspace for fft*/
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;

	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_real_workspace_alloc (n);
	real = gsl_fft_real_wavetable_alloc (n);

	/*Do the fft*/
	gsl_fft_real_transform (data_cp, 1, n, real, work);	
	
	/*Unpack the output into array with alternating real and imaginary part*/	
	gsl_fft_halfcomplex_unpack (data_cp, complex_coefficient,1,n);
	
	/*fill the output fft_data with the powerspectrum */
	for (i = 0; i < n; i++)
	{
		fft_data[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1])/n; 
	}
	
	/*Free memory of wavetable and workspace*/
	gsl_fft_real_wavetable_free (real);
	gsl_fft_real_workspace_free (work);
}


/* Shifts the input fft_data to center the 0 frequency */
void fft_shift(double *fft_data, int n) /* input data, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	
	/* make copy of fft_data as reference for the shift */ 
	double *fft_cp = malloc(n * sizeof (double));
	for (i = 0; i < n; i++)
	{
		fft_cp[i] = fft_data[i];
	}

	/* make shift */
	for (i = 0; i < n; i++)
	{
		if (n % 2) /*if n odd*/
		{ 
			if (i<=(n-2)/2)
			{
				fft_data[i] = fft_cp[(i+(n+1)/2)];
			}
			else
			{
				fft_data[i] = fft_cp[(i+(n+1)/2)%(n)];
			}			
		}
		else
		{
			if (i<n/2)
			{
				fft_data[i] = fft_cp[i+n/2];
			}
			else
			{
				fft_data[i] = fft_cp[(i+n/2)%(n)];
			}			
		}
	}
}

/* Makes a frequency array fft_freq with frequency interval 1/(dt*n) */
void fft_freq(double *fft_freq, double dt, int n) /* output frequency array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output array with frequencies */
	for (i = 0; i < n; i++)
	{
		fft_freq[i] = i/dt/n;
	}
}

/* Makes a frequency array fft_freq with frequency interval 1/(dt*n) with a centered O frequency */
void fft_freq_shift(double *fft_freq, double dt, int n) /* output frequency array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output aaray with shifted frequencies */
	for (i = 0; i < n; i++)
	{
		if (n % 2) /*if n odd*/
		{ 
			if (i<=(n-2)/2)
			{
				fft_freq[i] = (-(n-1)/2+i)/dt/n;
			}
			else
			{
				fft_freq[i] = (i-(n-1)/2)/dt/n;
			}			
		}
		else
		{
			if (i<n/2)
			{
				fft_freq[i] = (-n/2+i)/dt/n;
			}
			else
			{
				fft_freq[i] = (i-n/2)/dt/n;
			}			
		}
	}
}


