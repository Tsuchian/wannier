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
	out = fopen("topowan.txt", "w");
	double kx, ky;
	int i, j;
	double M = 1.0;
	double  rx, ry;
	double wan[101][101][2] = { 0.0 };
	for (int x = 0; x < 101; x++)
	{
		for (int y = 0; y < 101; y++)
		{
			rx = -5.0 + double(x) / 10.0;
			ry = -5.0 + double(y) / 10.0;
			complex <double> v1={0.0,0.0},v2={0.0,0.0},v3={0.0,0.0},v4={0.0,0.0};
			for (int n = 0; n < 101; n++)
			{
				for (int m = 0; m < 101; m++)
				{
					kx = -pi + 2.0*pi* double(n) / 100.0;
					ky = -pi + 2.0*pi* double(m) / 100.0;
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
					double kr = kx*rx + ky*ry;
					complex <double> v[4] = { {0,0} };
					for (j = 0; j < 4; j++)
					{
						v[j] = { real(b[j])*cos(kr) + imag(b[j])*sin(kr),imag(b[j])*cos(kr) - real(b[j])*sin(kr) };
					}
					v1 = { real(v1) + real(v[0]) , imag(v1) + imag(v[0] };
					v2 = { real(v2) + real(v[1]) , imag(v2) + imag(v[1] };
					v3 = { real(v3) + real(v[2]) , imag(v3) + imag(v[2] };
					v4 = { real(v4) + real(v[3]) , imag(v4) + imag(v[3] };
					//v1=v1 + v1;
					//v1=v1 + v1;
					delete[]a;
				}
			}
			wan[x][y][0]=real(v1)*real(v1)+imag(v1)*imag(v1)+real(v2)*real(v2)+imag(v2)*imag(v2);
			wan[x][y][1]=real(v3)*real(v3)+imag(v3)*imag(v3)+real(v4)*real(v4)+imag(v4)*imag(v4);
			fprintf(out, "%f %f ", rx, ry);
			for (i = 0; i < 2 * N; i++)
			{
				fprintf(out, " %f ", wan[x][y][i]);
			}
			fprintf(out, "\n");//write the eigenvalue
		}
	}

}
