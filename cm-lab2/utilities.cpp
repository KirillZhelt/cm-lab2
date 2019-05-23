#include <iostream>
#include <fstream>
#include <string>

#include "utilities.h"
#include "qr.h"

using namespace std;

ostream& operator<<(ostream& out, Complex complex_number) {
	out << complex_number.real;
	
	if (complex_number.imaginary != 0)
		out << " + " << complex_number.imaginary << "i";

	return out;
}

double CheckEigenvalue(double** A, int rows, int columns, double eigenvalue, double* eigenvector) {
	double* difference = new double[columns];
	double* Au = new double[columns];

	Multiply(A, rows, columns, eigenvector, columns, Au);

	for (int i = 0; i < columns; i++)
		difference[i] = Au[i] - eigenvalue * eigenvector[i];

	double euclidean_norm = EuclideanNorm(difference, columns);

	delete[] Au;
	delete[] difference;
	
	return euclidean_norm;
}

double MaxNorm(double* v, int length) {
	double max = v[0];

	for (int i = 1; i < length; i++) {
		if (fabs(v[i]) > max)
			max = fabs(v[i]);
	}

	return max;
}

double EuclideanNorm(double* v, int length) {
	return sqrt(ScalarMultiply(v, v, length));
}

double ScalarMultiply(double* v1, double* v2, int length) {
	double result = 0;

	for (int i = 0; i < length; i++)
		result += v1[i] * v2[i];

	return result;
}

void PrintMatrix(double** m, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			cout << m[i][j] << ' ';

		cout << '\n';
	}
}

void WriteMatrixToFile(double** m, int rows, int columns, string filename) {
	ofstream fout(filename);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			fout << m[i][j] << " ";
		}

		fout << "\n";
	}

	fout.close();
}

void CopyMatrix(double** src, double** dst, int rows, int columns) {
	for (int i = 0; i < rows; i++)
		memcpy(dst[i], src[i], columns * sizeof(double));
}

void CopyVector(double* src, double* dst, int length) {
	memcpy(dst, src, length * sizeof(double));
}

void Multiply(double** m, int rows, int columns,
	double* v, int length, double* b) {

	for (int i = 0; i < rows; i++) {
		b[i] = 0;

		for (int j = 0; j < columns; j++)
			b[i] += m[i][j] * v[j];
	}
}

void MatrixMultiply(double** a, double** b, int n, int m, int k, double** c) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < k; j++) {
			c[i][j] = 0;

			for (int q = 0; q < m; q++) {
				c[i][j] += a[i][q] * b[q][j];
			}
		}
	}
}

void Subtract(double* v1, double* v2, int length, double* result) {
	for (int i = 0; i < length; i++)
		result[i] = v1[i] - v2[i];
}