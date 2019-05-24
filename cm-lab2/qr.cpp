#include <algorithm>

#include "utilities.h"

#include "qr.h"

using namespace std;

pair<Complex, Complex> SolveQuadraticEquation(double a, double b, double c) {
	pair<Complex, Complex> roots;

	double discriminant = b * b - 4 * a * c;

	if (discriminant < 0) {
		roots.first.real = -b / (2 * a);
		roots.first.imaginary = sqrt(-discriminant) / (2 * a);

		roots.second.real = -b / (2 * a);
		roots.second.imaginary = -sqrt(-discriminant) / (2 * a);
	}
	else if (discriminant == 0) {
		roots.first.real = (-b + sqrt(discriminant)) / (2 * a);
		roots.first.imaginary = 0;

		roots.second.real = (-b - sqrt(discriminant)) / (2 * a);
		roots.second.imaginary = 0;
	}
	else {
		// discriminant > 0
		roots.first.real = (-b + sqrt(discriminant)) / (2 * a);
		roots.first.imaginary = 0;

		roots.second.real = (-b - sqrt(discriminant)) / (2 * a);
		roots.second.imaginary = 0;
	}

	return roots;
}

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

		if (norm != 0) {
			for (int j = i; j < rows; j++) // wi
				qr[i][j] /= norm;

			for (int j = i + 1; j < columns; j++) {
				double scalar_multiply_result = ScalarMultiply(&qr[j][i], &qr[i][i], rows - i);

				for (int k = i; k < rows; k++)
					qr[j][k] -= 2 * scalar_multiply_result * qr[i][k];
			}
		}
	}
}

void BuildHessenbergArnoldi(double** A, int rows, int columns, double** H) {
	double** Q = new double*[rows];
	for (int i = 0; i < rows; i++)
		Q[i] = new double[columns];

	double* z = new double[rows];
	double* hq = new double[rows];

	for (int i = 0; i < rows; i++) {
		memset(H[i], 0, sizeof(double) * columns);
	}

	double euclidean_norm = 0;
	for (int i = 0; i < rows; i++) {
		euclidean_norm += A[i][0] * A[i][0];
	}

	euclidean_norm = sqrt(euclidean_norm);

	for (int i = 0; i < rows; i++) {
		Q[i][0] = A[i][0] / euclidean_norm; // ||q1|| = 1
	}

	for (int k = 1;; k++) {
		for (int i = 0; i < rows; i++) {
			z[i] = 0;

			for (int j = 0; j < columns; j++)
				z[i] += A[i][j] * Q[j][k - 1];
		}

		for (int i = 0; i < k; i++) {
			H[i][k - 1] = 0;

			for (int j = 0; j < rows; j++)
				H[i][k - 1] += z[j] * Q[j][i];

			for (int j = 0; j < rows; j++)
				hq[j] = H[i][k - 1] * Q[j][i];

			Subtract(z, hq, rows, z);
		}

		PrintMatrix(H, rows, columns);
		cout << endl;
		double l = EuclideanNorm(z, rows);
		H[k][k - 1] = l;
		if (H[k][k - 1] < 1e-13)
			break;

		for (int j = 0; j < rows; j++)
			Q[j][k] = z[j] / H[k][k - 1];
	}

	delete[] hq;

	delete[] z;

	for (int i = 0; i < rows; i++)
		delete[] Q[i];

	delete[] Q;
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

	double** ak = new double*[rows + 1];
	for (int i = 0; i < rows + 1; i++)
		ak[i] = new double[columns + 1];

	for (int i = 0; i < rows; i++)
		q[i][i] = 1;

	BuildHessenbergArnoldi(A, rows, columns, ak);
	PrintMatrix(ak, rows, columns);

	for (int i = 0; i < number_of_iterations; i++) {
		BuildQR(ak, rows, columns, qr, diag_r);

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

		MatrixMultiply(r, q, rows, columns, columns, ak);

		for (int i = 0; i < rows; i++) {
			memset(q[i], 0, sizeof(double) * columns);
			memset(r[i], 0, sizeof(double) * columns);
		}

		for (int i = 0; i < rows; i++)
			q[i][i] = 1;

		for (int i = 0; i < rows; i++) {
			for (int j = i; j < columns; j++) {
				if (abs(ak[j][i]) < ZERO_EPS)
					ak[j][i] = 0;
			}
		}
	}
	
	for (int i = 0; i < rows; ) {
		if (i + 1 < rows) {
			if (ak[i + 1][i] != 0) {
				auto roots = SolveQuadraticEquation(1, -ak[i + 1][i + 1] - ak[i][i],
					ak[i][i] * ak[i + 1][i + 1] - ak[i][i + 1] * ak[i + 1][i]);
				eigenvalues[i] = roots.first;
				eigenvalues[i + 1] = roots.second;

				i += 2;
			}
			else {
				eigenvalues[i] = { ak[i][i], 0 };
				i++;
			}
		} 
		else {
			eigenvalues[i] = { ak[i][i], 0 };
			i++;
		}
	}

	for (int i = 0; i < rows; i++)
		delete[] ak[i];

	delete[] ak;

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