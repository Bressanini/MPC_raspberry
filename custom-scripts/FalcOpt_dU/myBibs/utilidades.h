#include <stdio.h>		// printf(), fopen(), fprintf(), fscanf();
#include <stdlib.h>		// malloc();
#include <math.h> 		// ceil()
#include <time.h>		//
#include <cblas.h>	// cblas_dgemv(), CblasRowMajor, CblasNoTrans
//#include "utilidades.c"

void mat2csv(double* M, int NL, int NC, char NOME[20]);

double* mZero(int NL, int NC);

void AXmBU(double* A, double* X, double* B, double* U, int NX, int NU, int NIT);

void inicializa(double** A, double** B,  int* NX, int* NU, double* TS, int* NIT, double** U);

void csv2mat(double** M, int NL, int NC, char NOME[20]);
