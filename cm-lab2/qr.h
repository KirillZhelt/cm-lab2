#pragma once

struct Complex {
	double real;
	double imaginary;
};

void BuildQR(double** m, int rows, int columns, double** qr, double* diag_r);

void FindEigenvaluesQR(double** A, int rows, int columns, Complex* eigenvalues, int number_of_iterations);