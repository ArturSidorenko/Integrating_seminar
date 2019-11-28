#include "integrals.h"


double quad::simple_quad(REAL_FUNC f, const qu_formula &q, double a, double b) {
	int n = q.n;

    double ans = 0;
    double x;
    for(int i = 0; i < n; i++) {
        x = 0.5 * (a+b) + 0.5*(b-a) * q.points[i];
        ans += 0.5*(b-a)*q.weights[i] * f(x);
    }
    return ans;
}

double quad::comp_quad(REAL_FUNC f, const qu_formula &q, double a, double b, int N) {
    double ans = 0;
    double h = (b-a)/N;
    for(int i = 0; i < N; i++) {
        ans += quad::simple_quad(f, q, a + i * h, a + (i+1)*h);
    }
    return ans;
}

quad::qu_formula quad::Newton_Kotes() {
    int n = 4;
    double *w = new double[n], *p = new double [n];
    p[0] = -1;
	p[1] = -1. / 3;
	p[2] = 1. / 3;
	p[3] = 1;
   
	quad::find_weights(w, p); //calculating appropriate weights
	
    return qu_formula(n,w,p);

}

quad::qu_formula quad::Gauss() {
	int n = 4;
	double *w = new double[n], *p = new double[n];

	double rad1 = (15 + 2 * sqrt(30))/35, rad2 = (15 - 2 * sqrt(30))/35;
	
	p[0] = -sqrt(rad1);
	p[1] = -sqrt(rad2);
	p[2] =  sqrt(rad2);
	p[3] =  sqrt(rad1);

	quad::find_weights(w, p); //calculating appropriate weights

	return qu_formula(n, w, p);

}

quad::qu_formula quad::Naive() {
	int n = 1;
	double *w = new double[n], *p = new double[n];
	w[0] = 2;
	p[0] = 0;

	return qu_formula(n, w, p);

}

void quad::Lin_eq(double *ans, const double *A, const double *b){
	double det = A[0] * A[3] - A[1] * A[2];
	if (det == 0) throw std::invalid_argument("The matrix is singular");
	ans[0] = (b[0] * A[3] - A[1] * b[1]) / det;
	ans[1] = (A[0] * b[1] - b[0] * A[2]) / det;

}
void quad::find_weights(double * ans, const double * p) //weights finder. works only for 4 symmetric points
{
	int n = 2;
	double *A = new double[n*n];
	double *b = new double[n];
	double *t = new double[n];
	
	for (int j = 0; j < n; j++) A[j] = 1;

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) A[n*i+j] = A[n*(i-1) + j]*p[j]*p[j]; //Vanderwonde matrix
	}

	for (int i = 0; i < n; i++) {
		b[i] = (2. / (2*i + 1)); //the integral of x^2i over [-1,1]
	}

	std::cout << "Dingo\n";
	Lin_eq(t, A, b);
	ans[0] = ans[3] = t[0]*0.5;
	ans[1] = ans[2] = t[1]*0.5;

	std::cout << "hello\n";

	delete[] A;
	delete[] b;
	delete[] t;
}



