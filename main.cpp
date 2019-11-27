#include "integrals.h"



using namespace std;
using namespace quad;

double func(double x) {return x*x*x*x;}

int main() {
	qu_formula q = Newton_Kotes(1);
    cout << comp_quad(*func, q, -1, 1, 346)<<"\n";
	system("pause");
    return 0;
}
