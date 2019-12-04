#include "integrals.h"

using namespace std;
using namespace quad;
using namespace qu2d;

double func(double x) { return x*x*x*x*x*x*x; }
double func2d(double x, double y) { return x*x*x*x*x*x*x; }

void test_lin_sys();
void test_qu_gen();
void test_2d_int();
void test_2d_partition(); 

void integral_comparison_1d();
void integral_comparison_2d();

const double PI = atan(1) * 4;
double lin_func(double x, double y) {
	return x + y;
}

int main() {
    //test_lin_sys();
	test_qu_gen();
    qu_formula q = Gauss();
    //cout << "The integral of x^7 over [0,1] equals "<<comp_quad(*func, q, 0, 1, 1)<<"\n";

	qu2d::partition p(1, 1, qu2d::point(0, 0), qu2d::point(1, 1));
	cout << "The integral of x+y  equals " << p.integrate(lin_func) << "\n";

    //test_2d_int();
	//integral_comparison_1d();
	integral_comparison_2d();

	//test_2d_partition();

	system("pause");
    return 0;
}

void test_lin_sys() {
	double *A = new double[4];
	double *b = new double[2];

        A[0] = 1; A[1] = 0; A[2] = 0.5; A[3] = 1;
	b[0] = 10; b[1] = 15;


	try {
		quad::Lin_eq(2, A, b);
                std::cout << "\n" << "Solutions: " << b[0] << " " << b[1] << "\n";
	}
	catch (std::invalid_argument &a) {
                std::cerr << a.what()<<"\n";
	}

	delete[] A;
	delete[] b;
}

void test_qu_gen() {
	qu_formula q (Naive()), b;
	qu_formula a(q);
	q = b = a;
	q = q;
	if (fabs(comp_quad(func, q, 0, 1, 100) - comp_quad(func, b, 0, 1, 100)) > 0.001) std::cout << "oups\n";
}

void test_2d_partition() {
	point a(0, 0), b(1, 1);
	qu2d::partition p, q, r;
	p = qu2d::partition(200, 200, a, b);
	q = qu2d::partition(200, 2000, a, b);
	r = p = q;
	if (fabs(r.integrate(func2d) - q.integrate(func2d)) > 0.001) std::cout << "oups\n";
}

void test_2d_int() {
	point a(0, 0), b(1, 1);
	qu2d::partition p(1000, 1000, a, b);

	p.save("test.txt");

	qu2d::partition q("test.txt");

	cout <<"The integral over the unit square before saving: " <<p.integrate(func2d)<<"\n";
	cout << "The integral over the unit square after saving: " << q.integrate(func2d) << "\n";
}

double exp_func(double x) {
	return exp(x) * (x * x - 2 * x + 0.5);
}

double highly_osc_func(double x) {
	return sin(100.27 * PI * x);
}

double discont_func(double x) {
	return (x < 1e-16) ? 0 : 1./sqrt(x);
}

void integral_comparison_1d() {
	qu_formula naive(Naive()), nc(Newton_Kotes()), g(Gauss());
	ofstream out("int_1d_pol.csv");

	const double func_val = 0.125, exp_val = -0.4225772573114326, osc_val = 0.0010751748439276885,
				 disc_val = 99999999. / 50000000;
	const int N = 10000;

	//pol func
	out << "n, naive, nc, gauss\n";
	for (int n = 1; n < N; n *= 2) {
		out << n << " , " << fabs(comp_quad(func, naive, 0, 1, n) - func_val) << " , " <<
			fabs(comp_quad(func, nc, 0, 1, n) - func_val) << " , " <<
			fabs(comp_quad(func, g, 0, 1, n) - func_val) << "\n";
	}

	out.close();
	out = ofstream("int_1d_exp.csv");

	//exp func
	out << "n, naive, nc, gauss\n";
	for (int n = 1; n < N; n *= 2) {
		out << n << " , " << fabs(comp_quad(exp_func, naive, 0, 1, n) - exp_val) << " , " <<
			fabs(comp_quad(exp_func, nc, 0, 1, n) - exp_val) << " , " <<
			fabs(comp_quad(exp_func, g, 0, 1, n) - exp_val) << "\n";
	}
	out.close();

	std::cout << "Exp int = " << comp_quad(exp_func, g, 0, 1, N) << "\n";

	out = ofstream("int_1d_sin.csv");
	//sine func
	out << "n, naive, nc, gauss\n";
	for (int n = 1; n < N; n *= 2) {
		out << n << " , " << fabs(comp_quad(highly_osc_func, naive, 0, 1, n) - osc_val) << " , " <<
			fabs(comp_quad(highly_osc_func, nc, 0, 1, n) - osc_val) << " , " <<
			fabs(comp_quad(highly_osc_func, g, 0, 1, n) - osc_val) << "\n";
	}

	out.close();
	out = ofstream("int_1d_sqrt.csv");
	//discontinious function
	out << "n, naive, nc, gauss\n";
	for (int n = 1; n < N; n *= 2) {
		out << n << " , " << fabs(comp_quad(discont_func, naive, 0, 1, n) - disc_val) << " , " <<
			fabs(comp_quad(discont_func, nc, 0, 1, n) - disc_val) << " , " <<
			fabs(comp_quad(discont_func, g, 0, 1, n) - disc_val) << "\n";
	}

	out.close();

}



double exp_func(double x, double y) {
	return exp(x + 2*y) * (2*x - y);
}

double sin_func(double x, double y) {
	return sin(100.27 * PI * (x - 2*y));
}

double disc_func(double x, double y) {
	return 1./(sqrt(x) + sqrt(y));
}


void integral_comparison_2d() {
	ofstream out("int_2d.csv");
	qu2d::point a(0, 0), b(1, 1);
	;

	double exp_val = 2.785365435751634;
	double sin_val = -2.5602582924483034e-6;
	double disc_val = 0.8182741851734792;
	

	out << "n,exp,sin\n";
	for (int n = 1; n < 1000; n*=2) {
		qu2d::partition p(n, n, a, b);
		out << n << "," << fabs(p.integrate(exp_func) - exp_val) << "," << fabs(p.integrate(sin_func) - sin_val) <<  "\n";
	}
	out.close();

}