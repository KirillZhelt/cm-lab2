#pragma once

struct Complex {
	double real;
	double imaginary;
};

void FindEigenvaluesQR(double** A, int rows, int columns, Complex* eigenvalues);