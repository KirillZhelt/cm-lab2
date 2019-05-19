
#include <iostream>

#include "utilities.h"

#include "power_iteration.h"

using namespace std;

void PowerIteration(double** A, int rows, int columns, double& eigenvalue1, double& eigenvalue2,
	double* eigenvector1, double* eigenvector2) {
	double* y = new double[columns] {};
	y[0] = 1;

	eigenvalue1 = 1;

	CopyVector(y, eigenvector1, columns);

	while (true) {

	}

	delete[] y;
}