#include "utilities.h"

#include "qr.h"

using namespace std;

int Sign(double d) {
	return (d >= 0) - (d < 0);
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

void FindEigenvaluesQR(double** A, int rows, int columns, Complex* eigenvalues, int number_of_iterations) {
	// Ak = QkRk
	// Ak+1 = Rk*Qk
	// Считаю собственные значения из Ak

	double** qr = new double*[rows];
	for (int i = 0; i < rows; i++)
		qr[i] = new double[columns];

	double* diag_r = new double[rows];

	double** r = new double*[rows];
	for (int i = 0; i < rows; i++)
		r[i] = new double[columns] {};

	double** q = new double*[rows];
	for (int i = 0; i < rows; i++)
		q[i] = new double[columns] {};

	for (int i = 0; i < rows; i++)
		q[i][i] = 1;

	for (int i = 0; i < number_of_iterations; i++) {
		BuildQR(A, rows, columns, qr, diag_r);

		for (int i = 0; i < rows; i++)
			r[i][i] = diag_r[i];

		for (int i = 0; i < rows; i++) {
			for (int j = i + 1; j < columns; j++)
				r[i][j] = qr[j][i];
		}

		for (int i = 0; i < rows; i++) {
			// every row go through w1 ... wn

			for (int j = 0; j < columns - 1; j++) {
				double scalar_multiply_result = ScalarMultiply(&qr[j][j], &q[i][j], rows - j);
				
				for (int k = j; k < columns; k++)
					q[i][k] -= 2 * scalar_multiply_result * qr[j][k];
			}
		}
	}

	/*
		https://www.programiz.com/cpp-programming/examples/quadratic-roots

		for finding roots 
	*/

	for (int i = 0; i < rows; i++)
		delete[] q[i];

	delete[] q;

	for (int i = 0; i < rows; i++)
		delete[] r[i];

	delete[] r;

	for (int i = 0; i < rows; i++)
		delete[] qr[i];

	delete[] qr;
}