#pragma once

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<exception>
#include<string>

namespace quad {

// this thing stores all th coefficients
//class structure prevents any changes during calculations
//
//NB this class stores data for the segment [-1, 1]

class LogicError
{
private:
    std::string s;
public:
	LogicError(std::string r) { s = r; };
	std::string get() { return s; };
};

struct qu_formula
{
    int n;
    double *weights, *points;
	
	qu_formula() {};
	qu_formula(int n_, double *w_, double *p_) : n(n_), weights(w_), points(p_) {};
};

typedef  double (*REAL_FUNC)(double); //the real function

double simple_quad(REAL_FUNC f, const qu_formula &q, double a, double b); //simple quadrature formula
double comp_quad(REAL_FUNC f, const qu_formula &q, double a, double b, int N); //compound quadrature formula

qu_formula Newton_Kotes(int n); //making quadratures
qu_formula Gauss(int n);

void Lin_eq(double *ans, const double *A, const double *b); //Linear system solving

}

