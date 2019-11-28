#pragma once

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<exception>
#include<string>
#include<algorithm>
#include<utility>
#include<cstring>

namespace quad {

	// this thing stores all th coefficients
	//class structure prevents any changes during calculations
	//
	//NB this class stores data for the segment [-1, 1]


	class qu_formula
	{
    private:
		int n_;
		double *weights_, *points_;
    public:
		qu_formula();
		qu_formula(int n, double *w, double *p) : n_(n), weights_(w), points_(p) {};
        qu_formula(const qu_formula &q); //copy constructor
        qu_formula(qu_formula &&q); //move constructor
        const qu_formula& operator=(const qu_formula &q);
        const qu_formula& operator=(qu_formula &&other);
        ~qu_formula();
        
        int get_n() {return n_;};
        double get_weight(int k);
        double get_point(int k);
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
