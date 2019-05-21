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
		for (int j = 0; j < columns; j++)
			m[i][j] = dis(gen);
	}
}
