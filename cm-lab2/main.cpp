
#include <iostream>
#include <fstream>
#include <chrono>

#include "utilities.h"

#include "fill.h" // TASK1
#include "power_iteration.h" // TASK2
#include "qr.h" // TASK3

using namespace std;

const int N = 3;

const int ROWS = 10;
const int COLUMNS = 10;

const int NUMBER_OF_ITERATIONS_QR = 1000;

int main() {
	ofstream fout("report.txt");

	chrono::high_resolution_clock::time_point start, finish;
	chrono::duration<double, std::milli> fp_ms;

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
	/*
	start = chrono::high_resolution_clock::now();
	PowerIteration(A, ROWS, COLUMNS, eigenvalue1, eigenvalue2, eigenvector1, eigenvector2, k1, k2);
	finish = chrono::high_resolution_clock::now();

	fp_ms = finish - start;

	fout << "POWER ITERATION (TASK 2):" << endl;
	fout << "Average time to count one eigenvalue: " << fp_ms.count() / 2 << " ms" << endl;
	fout << "Number of iterations for first and second eigenvalues: " << k1 << " " << k2 << endl;
	fout << "First eigenvalue norm: " << CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue1, eigenvector1) << endl;
	fout << "Second eigenvalue norm: " << CheckEigenvalue(A, ROWS, COLUMNS, eigenvalue2, eigenvector2) << endl;
	fout << endl << endl;
	*/
	// QR (TASK 3)
	Complex* eigenvalues = new Complex[ROWS];

	start = chrono::high_resolution_clock::now();
	FindEigenvaluesQR(A, ROWS, COLUMNS, eigenvalues, NUMBER_OF_ITERATIONS_QR);
	finish = chrono::high_resolution_clock::now();

	fp_ms = finish - start;

	fout << "QR (TASK 3): " << endl;
	fout << "Average time to count all eigenvalues: " << fp_ms.count() << " ms" << endl;
	fout << "Eigenvalues: ";
	PrintVector(fout, eigenvalues, ROWS);
	fout << endl << "Numpy eigenvalues: " << endl;
	system("python eigenvalues.py"); // будет работать только с питоном и нумпаем
	/*
	// TASKS 5-8
	system("python equation.py"); // будет работать только с питоном
	*/
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