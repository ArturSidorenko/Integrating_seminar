#include 'integrals.h'


quad::qu_formula(int n) {
    n_ = n;
    weights_ = new double[n];
    points_ = new double[n];
}

quad::qu_formula(int n, double *w, double *p) {
    n_ = n;
    weights = w_;
    points = p_;
}

quad::~qu_formula() {
    delete [] points;
    delete [] weights;
}

double quad::get_point(int i) {
    if((i < 0) || (i > n_)) throw quad::LogicError("Invalid massive index");
    return points[i];
}

double quad::get_weight(int i) {
    if((i < 0) || (i > n_)) throw quad::LogicError("Invalid massive index");
    return weights[i];
}

double quad::simple_quad(REAL_FUNC f, const &qu_formula q, double a, double b) {
    int n = q.get_n();

    double ans = 0;
    double x;
    for(int i = 0; i < n; i++) {
        x = 0.5 * (a+b) + 0.5*(b-a) * q.get_point(i);
        ans += 0.5*(b-a)*q.get_weight(i) * f(x);
    }
    return ans;
}

double quad::comp_quad(REAL_FUNC f, const &qu_formula q, double a, double b, int N) {
    double ans = 0;
    double h = (b-a)/N;
    for(int i = 0; i < N; i++) {
        ans += quad::simple_quad(f, q, a + i * h, a + (i+1)*h);
    }
    return ans;
}

qu_formula quad::Newton_Kotes(int n) {
    n = 1;
    double *w = new double[n], *p = new double [n];
    w[0] = 2;
    p[0] = 0;

    return qu_formula(n,w,p);

}

qu_formula quad::Gauss(int n) {
    n = 1;
    double *w = new double[n], *p = new double [n];
    w[0] = 2;
    p[0] = 0;

    return qu_formula(n,w,p);

}

void quad::Lin_eq(double *ans, const double *A, const double *b){};


