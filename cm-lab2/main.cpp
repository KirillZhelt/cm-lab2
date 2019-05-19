
#include <iostream>

#include "utilities.h"

#include "fill.h" // TASK1
#include "power_iteration.h" // TASK2

using namespace std;

const int N = 3;

const int ROWS = 10;
const int COLUMNS = 10;

int main() {
	double** A = new double*[ROWS];

	for (int i = 0; i < ROWS; i++)
		A[i] = new double[COLUMNS];

	Fill(A, ROWS, COLUMNS, N);
	PrintMatrix(A, ROWS, COLUMNS);

	double eigenvalue1, eigenvalue2;

	double* eigenvector1 = new double[COLUMNS];
	double* eigenvector2 = new double[COLUMNS];

	PowerIteration(A, ROWS, COLUMNS, eigenvalue1, eigenvalue2, eigenvector1, eigenvector2);

	cout << endl << endl;
	cout << CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue1, eigenvector1) << endl;
	cout << CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue2, eigenvector2);

	delete[] eigenvector2;
	delete[] eigenvector1;

	for (int i = 0; i < ROWS; i++)
		delete[] A[i];

	delete[] A;

	system("pause");

	return 0;
}