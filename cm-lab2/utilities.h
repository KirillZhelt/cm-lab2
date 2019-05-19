#pragma once

double MaxNorm(double* v, int length);

double EuclideanNorm(double* v, int length);

double ScalarMultiply(double* v1, double* v2, int length);

void PrintMatrix(double** m, int rows, int columns);

void PrintVector(double* v, int length);

void CopyMatrix(double** src, double** dst, int rows, int columns);

void CopyVector(double* src, double* dst, int length);

void Multiply(double** m, int rows, int columns,
	double* v, int length, double* b);

void Subtract(double* v1, double* v2, int length, double* result);

