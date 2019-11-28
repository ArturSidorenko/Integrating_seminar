#pragma once

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<exception>
#include<string>
#include<algorithm>
#include<utility>

namespace quad {

	// this thing stores all th coefficients
	//class structure prevents any changes during calculations
	//
	//NB this class stores data for the segment [-1, 1]


	struct qu_formula
	{
		int n;
		double *weights, *points;

		qu_formula() {};
		qu_formula(int n_, double *w_, double *p_) : n(n_), weights(w_), points(p_) {};
	};

	typedef  double(*REAL_FUNC)(double); //the real function

	double simple_quad(REAL_FUNC f, const qu_formula &q, double a, double b); //simple quadrature formula
	double comp_quad(REAL_FUNC f, const qu_formula &q, double a, double b, int N); //compound quadrature formula

	qu_formula Newton_Kotes(); //making quadratures
	qu_formula Gauss();
	qu_formula Naive();

	double col_max(double *A, int n, int k); //to control errors
	void swap_rows(double *A, int n, int p, int q);
	void add_rows(double *A, int n, int p, int q, double m);
	void divide_by(double *A, int n, int p, double d);
	void Lin_eq(int n, double *A, double *b);  //Linear system solving. The matrix A will be spoiled, the answer will be in the vector b

	void find_weights(int n, double *ans, const double *p);

}