#pragma once

const double POWER_ITERATION_EPS = 0.00001;
const int MAX_ITERATION_NUMBER = 500;

void PowerIteration(double** A, int rows, int columns, double& eigenvalue1, double& eigenvalue2, 
	double* eigenvector1, double* eigenvector2);