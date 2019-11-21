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
    string s;
public:
    LogicError(string r) {s = r};
    string get() {return s;}
}

class qu_formula
{
private:
    int n_;
    double *weights_, *points_;
public:
    qu_formula(int n);
    qu_formula(int n, double *w, double *p);
    ~qu_formula();
    int get_n() {return n_;}
    double get_point(int i);
    double get_weight(int i);

};

typedef  double (*REAL_FUNC)(double); //the real function

double simple_quad(REAL_FUNC f, const &qu_formula q, double a, double b); //simple quadrature formula
double comp_quad(REAL_FUNC f, const &qu_formula q, double a, double b, int N); //compound quadrature formula

qu_formula Newton_Kotes(int n); //making quadratures
qu_formula Gauss(int n);

void Lin_eq(double *ans, const double *A, const double *b); //Linear system solving

}

