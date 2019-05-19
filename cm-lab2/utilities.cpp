#include <iostream>

#include "utilities.h"

using namespace std;

void PrintMatrix(double** m, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++)
			cout << m[i][j] << ' ';

		cout << '\n';
	}
}