#include <random>
#include <iostream>

#include "fill.h"

using namespace std;

void Fill(double** m, int rows, int columns, int N) {
	double limit = pow(2, N / 4.0);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dis(-limit, limit);

	for (int i = 0; i < rows; i++)
		m[i][i] = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = i + 1; j < columns; j++) {
			m[i][j] = m[j][i] = dis(gen);

			m[i][i] += fabs(m[i][j]);
			m[j][j] += fabs(m[i][j]);
		}

		if (m[i][i] < 0)
			cout << m[i][i] << endl;
	}
}
