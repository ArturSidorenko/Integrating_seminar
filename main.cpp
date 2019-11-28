#include "integrals.h"



using namespace std;
using namespace quad;

double func(double x) { return x*x*x*x*x*x*x; }

void test_lin_sys();

const extern double PI = atan(1) * 4;

int main() {
	test_lin_sys();
	qu_formula q = Gauss();
    cout << comp_quad(*func, q, 0, 1, 1)<<"\n";

	system("pause");
    return 0;
}

void test_lin_sys() {
	double *A = new double[4];
	double *b = new double[2];

	A[0] = 1; A[1] = 0; A[2] = 0; A[3] = 0;
	b[0] = 10; b[1] = 15;


	try {
		quad::Lin_eq(2, A, b);
	}
	catch (std::invalid_argument &a) {
		std::cerr << a.what();
	}
	std::cout << "\n" << "Solutions: " << b[0] << " " << b[1] << "\n";
	delete[] A;
	delete[] b;
}