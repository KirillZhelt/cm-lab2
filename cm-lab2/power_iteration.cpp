
#include <iostream>

#include "utilities.h"

#include "power_iteration.h"

using namespace std;

void PowerIterationBase(double** A, int rows, int columns, double& eigenvalue, double* eigenvector, int& k) {
	double* y = new double[columns] {};
	y[0] = 1;

	double* difference = new double[columns];

	double* Au = new double[columns];

	eigenvalue = 0;

	CopyVector(y, eigenvector, columns);

	k = 0;

	Multiply(A, rows, columns, eigenvector, columns, Au);

	while (k < MAX_ITERATION_NUMBER) {

		for (int i = 0; i < columns; i++)
			difference[i] = Au[i] - eigenvalue * eigenvector[i];

		if (EuclideanNorm(difference, columns) < POWER_ITERATION_EPS)
			break;

		CopyVector(Au, y, columns);

		double y_norm = EuclideanNorm(y, columns);

		for (int i = 0; i < columns; i++)
			eigenvector[i] = y[i] / y_norm;

		Multiply(A, rows, columns, eigenvector, columns, Au);

		eigenvalue = ScalarMultiply(eigenvector, Au, columns);

		k++;
	}

	delete[] Au;

	delete[] difference;

	delete[] y;
}

void PowerIteration(double** A, int rows, int columns, double& eigenvalue1, double& eigenvalue2,
	double* eigenvector1, double* eigenvector2, int& k1, int& k2) {
	PowerIterationBase(A, rows, columns, eigenvalue1, eigenvector1, k1);

	double** B = new double*[rows];
	for (int i = 0; i < rows; i++)
		B[i] = new double[columns];

	CopyMatrix(A, B, rows, columns);

	for (int i = 0; i < rows; i++)
		B[i][i] -= eigenvalue1;

	PowerIterationBase(B, rows, columns, eigenvalue2, eigenvector2, k2);

	eigenvalue2 += eigenvalue1;

	for (int i = 0; i < rows; i++)
		delete[] B[i];

	delete[] B;
}