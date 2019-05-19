
#include <iostream>

#include "utilities.h"

#include "power_iteration.h"

using namespace std;

void PowerIteration(double** A, int rows, int columns, double& eigenvalue1, double& eigenvalue2,
	double* eigenvector1, double* eigenvector2) {
	double* y = new double[columns] {};
	y[0] = 1;

	double* difference = new double[columns];

	double* Au = new double[columns];

	eigenvalue1 = 0;

	CopyVector(y, eigenvector1, columns);

	int k = 0;
	
	Multiply(A, rows, columns, eigenvector1, columns, Au);

	while (k < MAX_ITERATION_NUMBER) {

		for (int i = 0; i < columns; i++)
			difference[i] = Au[i] - eigenvalue1 * eigenvector1[i];

		if (EuclideanNorm(difference, columns) < POWER_ITERATION_EPS)
			break;

		CopyVector(Au, y, columns);

		double y_norm = EuclideanNorm(y, columns);

		for (int i = 0; i < columns; i++)
			eigenvector1[i] = y[i] / y_norm;

		Multiply(A, rows, columns, eigenvector1, columns, Au);

		eigenvalue1 = ScalarMultiply(eigenvector1, Au, columns);

		k++;
	}

	cout << endl << k << endl;

	delete[] Au;

	delete[] difference;

	delete[] y;
}