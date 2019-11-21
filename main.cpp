#include 'integrals.h'

using namespace std;
using namespace quad;

double func(x) {return 1;}

int main() {
    qu_formula q = Newton_Kotes(1)
    cout << simple_quad(*func, q, -1, 1);
    return 0;
}
