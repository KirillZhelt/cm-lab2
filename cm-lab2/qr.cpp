#include "utilities.h"

#include "qr.h"

using namespace std;

int Sign(double d) {
	return (d >= 0) - (d < 0);
}

void Transpose(double** m, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = i + 1; j < columns; j++)
			swap(m[i][j], m[j][i]);
	}
}

void BuildQR(double** m, int rows, int columns, double** qr, double* diag_r) {
	// copy transpose

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			qr[i][j] = m[j][i];
	}

	for (int i = 0; i < columns; i++) {
		/*
			1) итерируемся по каждому столбцу
				1a) находим new_a = -signa1 * ||a|| (diag_r[i] = new_a[i][i])
				1b) находим w = (a - new_a) / ||a - new_a||
				1c) для каждого следующего столбца
					1c-a)  new_a[i] = a[i] - 2 * (a[i], w) * w
		*/

		diag_r[i] = Sign(qr[i][i]) * EuclideanNorm(&qr[i][i], rows - i);

		qr[i][i] -= diag_r[i];
		double norm = EuclideanNorm(&qr[i][i], rows - i);

		for (int j = i; j < rows; j++) // wi
			qr[i][j] /= norm;

		for (int j = i + 1; j < columns; j++) {
			double scalar_multiply_result = ScalarMultiply(&qr[j][i], &qr[i][i], rows - i);

			for (int k = i; k < rows; k++)
				qr[j][k] -= 2 * scalar_multiply_result * qr[i][k];
		}
	}
}

void SolveQR(double** qr, double* diag_r, int rows, int columns, double* v, double* x) {
	double* b = new double[rows];
	CopyVector(v, b, rows);

	for (int i = 0; i < columns - 1; i++) {
		double scalar_multiply_result = ScalarMultiply(&qr[i][i], &b[i], rows - i);

		for (int j = i; j < rows; j++)
			b[j] -= 2 * scalar_multiply_result * qr[i][j];
	}

	// Rx = b'
	for (int i = columns - 1; i >= 0; i--) {
		double sum = 0;

		for (int j = i + 1; j < columns; j++)
			sum += x[j] * qr[j][i];

		x[i] = (b[i] - sum) / diag_r[i];
	}


	delete[] b;
}

void FindEigenvaluesQR(double** A, int rows, int columns, Complex* eigenvalues) {

}