#include "integrals.h"



using namespace std;
using namespace quad;

double func(double x) { return x*x*x*x*x*x*x; }

void test_lin_sys();

const extern double PI = atan(1) * 4;

int main() {
	qu_formula q = Gauss();
    cout << comp_quad(*func, q, 0, 1, 1)<<"\n";

	system("pause");
    return 0;
}

void test_lin_sys() {
	double *A = new double[4];
	double *b = new double[2];
	double *ans = new double[2];

	A[0] = 1; A[1] = 1; A[2] = 2; A[3] = 1;
	b[0] = 10; b[1] = 15;

	quad::Lin_eq(ans, A, b);
	std::cout << "\n" << "Solutions: " << ans[0] << " " << ans[1] << "\n";
	delete[] A;
	delete[] b;
	delete[] ans;
}