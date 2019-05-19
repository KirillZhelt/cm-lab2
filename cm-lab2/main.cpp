
#include <iostream>

#include "utilities.h"

#include "fill.h" // TASK1

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

	for (int i = 0; i < ROWS; i++)
		delete[] A[i];

	delete[] A;

	system("pause");

	return 0;
}