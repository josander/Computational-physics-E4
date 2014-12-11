/*
fft_func.h
 
Created by Martin Gren on 2014-10-22.
*/

#ifndef _fft_func_h
#define _fft_func_h

extern void fft(double *, double *, int);

extern void fft_freq(double *, double, int);

extern void fft_shift(double *, int);

extern void fft_freq_shift(double *, double, int);

extern void make_fft_shift(double *, double *, double, int);


#endif
