#ifndef lapack_complex_double 
#define lapack_complex_double MKL_Complex16
#endif
#include <iostream>
#include <cstdlib>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include<fstream>
#include<cstring>
#define pi 3.1415926535898 
using namespace std;
int main()
{
	int N = 1;
	FILE *out, *matrix, *pm;
	out = fopen("topo.txt", "w");
	double kx, ky;
	int i, j;
	double M = 1.0;
	int t = 0;
	double k[200 + int(sqrt(2) * 100)][2];
	for (i = 0; i < 101; i++)
	{
		k[i][0] = pi*double(i) / 100;
		k[i][1] = 0.0;
	}
	for (i = 100; i < 201; i++)
	{
		k[i][0] = pi;
		k[i][1] = pi*double(i - 100) / 100;

	}
	for (i = 200; i<int(sqrt(2) * 100) + 200; i++)
	{
		k[i][0] = pi - pi*double(i - 200) / (sqrt(2) * 100);
		k[i][1] = pi - pi*double(i - 200) / (sqrt(2) * 100);
	}
	for (int n = 0; n < 200 + int(sqrt(2) * 100); n++)
	{
		kx = k[n][0];
		ky = k[n][1];
		//kx = -pi + 2.0*pi* double(n) / 100.0;
		//ky = -pi + 2.0*pi* double(m) / 100.0;
		complex<double> **a = new complex<double>*[2 * N];
		for (i = 0; i < 2 * N; i++)
		{
			a[i] = new complex<double>[2 * N];
		}
		a[0][0] = { M + 2.0 - cos(kx) - cos(ky),0.0 };
		a[0][1] = { sin(kx),-sin(ky) };
		a[1][0] = { sin(kx),sin(ky) };
		a[1][1] = { -(M + 2.0 - cos(kx) - cos(ky)),0.0 };
		complex<double> *b = new complex<double>[2 * N * 2 * N];
		for (int k = 0; k < 2 * N; k++)
		{
			for (int l = 0; l < 2 * N; l++)
			{
				b[k * 2 * N + l] = a[k][l];
			}
		}
		char jobz = 'V';
		char uplo = 'U';
		double w[2 * N];
		lapack_int lda = 2 * N;
		LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, 2 * N, b, lda, w); //diagonalization
		fprintf(out, "%d %f %f ", n, kx, ky);
		for (i = 0; i < 2 * N; i++)
		{
			//fprintf(out, "%f %f ",kx,ky);
			fprintf(out, " %f ", w[i]);
		}
		fprintf(out, "\n");//write the eigenvalue
		delete[]a;
	}
}
