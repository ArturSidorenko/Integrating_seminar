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

const double PI = atan(1) * 4;

int main() {
    //test_lin_sys();
	test_qu_gen();
    qu_formula q = Gauss();
    cout << "The integral of x^7 over [0,1] equals "<<comp_quad(*func, q, 0, 1, 1)<<"\n";


    test_2d_int();

	test_2d_partition();

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