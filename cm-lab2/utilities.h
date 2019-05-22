#pragma once

#include <iostream>

struct Complex;
std::ostream& operator<<(std::ostream& out, Complex complex_number);

double CheckEigenvalue(double** A, int rows, int columns, double eigenvalue, double* eigenvector);

double MaxNorm(double* v, int length);

double EuclideanNorm(double* v, int length);

double ScalarMultiply(double* v1, double* v2, int length);

void PrintMatrix(double** m, int rows, int columns);

void WriteMatrixToFile(double** m, int rows, int columns, std::string filename);

void CopyMatrix(double** src, double** dst, int rows, int columns);

void CopyVector(double* src, double* dst, int length);

void Multiply(double** m, int rows, int columns,
	double* v, int length, double* b);

void MatrixMultiply(double** a, double** b, int n, int m, int k, double** c);

void Subtract(double* v1, double* v2, int length, double* result);

template<typename T>
void PrintVector(T* v, int length) {
	for (int i = 0; i < length; i++)
		std::cout << v[i] <<  "     ";

	std::cout << "\n";
}
