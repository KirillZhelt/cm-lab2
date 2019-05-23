
#include <iostream>
#include <fstream>

#include "utilities.h"

#include "fill.h" // TASK1
#include "power_iteration.h" // TASK2
#include "qr.h" // TASK3

using namespace std;

const int N = 3;

const int ROWS = 10;
const int COLUMNS = 10;

int main() {
	ofstream fout("report.txt");

	double** A = new double*[ROWS];

	for (int i = 0; i < ROWS; i++)
		A[i] = new double[COLUMNS];

	// FILL (TASK 1)
	Fill(A, ROWS, COLUMNS, N);
	WriteMatrixToFile(A, ROWS, COLUMNS, "matrix.txt");

	// POWER ITERATION (TASK 2)
	double eigenvalue1, eigenvalue2;

	double* eigenvector1 = new double[COLUMNS];
	double* eigenvector2 = new double[COLUMNS];

	int k1, k2;

	PowerIteration(A, ROWS, COLUMNS, eigenvalue1, eigenvalue2, eigenvector1, eigenvector2, k1, k2);

	cout << k1 << " " << k2;

	cout << endl << endl;
	cout << eigenvalue1 << ": " << CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue1, eigenvector1) << endl;
	cout << eigenvalue2 << ": "<< CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue2, eigenvector2) << endl;
	cout << endl << endl;

	Complex* eigenvalues = new Complex[ROWS];

	FindEigenvaluesQR(A, ROWS, COLUMNS, eigenvalues, 1000);
	PrintVector(eigenvalues, ROWS);

	cout << endl << "Numpy eigenvalues: " << endl;
	system("python eigenvalues.py");
	cout << endl;


	delete[] eigenvalues;

	delete[] eigenvector2;
	delete[] eigenvector1;

	for (int i = 0; i < ROWS; i++)
		delete[] A[i];

	delete[] A;

	fout.close();

	system("pause");

	return 0;
}