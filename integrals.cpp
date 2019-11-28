#include "integrals.h"

quad::qu_formula::qu_formula() {
    n_ = 1;
    weights_ = new double[1];
    points_ = new double[1];
}

quad::qu_formula::qu_formula(const qu_formula &other) {
    n_ = other.n_;
    weights_ = new double[n_];
    points_ = new double[n_];
    std::memcpy(weights_, other.weights_, n_);
    std::memcpy(points_, other.points_, n_);
}

quad::qu_formula::~qu_formula() {
    delete [] points_;
    delete [] weights_;
}

const quad::qu_formula& quad::qu_formula::operator=(const qu_formula &other) {
    if(this != &other) {
        delete [] points_;
        delete [] weights_;
        n_ = other.n_;
        weights_ = new double[n_];
        points_ = new double[n_];
        std::memcpy(weights_, other.weights_, n_);
        std::memcpy(points_, other.points_, n_);
    }
    return *this;
}

quad::qu_formula::qu_formula(qu_formula &&other) {
    weights_ = other.weights_;
    points_ = other.points_;
    std::swap(n_, other.n_);
    other.weights_ = new double [1];
    other.points_ = new double [1];

}

const quad::qu_formula& quad::qu_formula::operator=(qu_formula &&other) {
    if(this !=&other) {
        delete [] points_;
        delete [] weights_;
        n_ = other.n_;
        weights_ = other.weights_;
        points_ = other.points_;
        other.weights_ = new double [1];
        other.points_ = new double [1];
    }
    return *this;
}

double quad::qu_formula::get_point(int k) {
    if(k>=n_) throw std::invalid_argument("The required point index is out of range");
    return points_[k];
}

double quad::qu_formula::get_weight(int k) {
    if(k>=n_) throw std::invalid_argument("The required weight index is out of range");
    return weights_[k];
}

double quad::simple_quad(REAL_FUNC f, const qu_formula &q, double a, double b) {
        int n = q.get_n();

    double ans = 0;
    double x;
    for(int i = 0; i < n; i++) {
        x = 0.5 * (a+b) + 0.5*(b-a) * q.get_point(i);
        ans += 0.5*(b-a)*q.get_weight(i) * f(x);
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
   
	quad::find_weights(n, w, p); //calculating appropriate weights
	
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

	quad::find_weights(n, w, p); //calculating appropriate weights

	return qu_formula(n, w, p);

}

quad::qu_formula quad::Naive() {
	int n = 1;
	double *w = new double[n], *p = new double[n];
	w[0] = 2;
	p[0] = 0;

	return qu_formula(n, w, p);

}


void quad::find_weights(int n, double * ans, const double * p) //weights finder. works only for 4 symmetric points
{
	double *A = new double[n*n];

	
	for (int j = 0; j < n; j++) A[j] = 1;

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) A[n*i+j] = A[n*(i-1) + j]*p[j]; //Vanderwonde matrix
	}

	for (int i = 0; i < n; i++) {
		ans[i] = (i%2==0) * (2. / (i + 1)); //the integral of x^i over [-1,1]
	}

	try {
		quad::Lin_eq(n, A, ans);
	}
	catch (const std::invalid_argument &ups) {
		std::cerr << "An invalid argument exception occured: "<< ups.what() << "\n";
		exit(-8);
	}
	catch (...) {
		std::cerr << "A critical error occured\n";
		exit(-1);
	}

	delete[] A;
}

void quad::swap_rows(double *A, int n, int p, int q) {
	for (int j = 0; j < n; j++) std::swap(A[n*p + j], A[n*q + j]);
}

void quad::add_rows(double *A, int n, int p, int q, double m) {
	for (int j = 0; j < n; j++) A[n*p + j] += m * A[n*q + j];
}
void quad::divide_by(double *A, int n, int p, double d) {
	for (int j = 0; j < n; j++) A[n*p + j] /= d;
}

double quad::col_max(double *A, int n, int k) {
	double m = 0;
	for (int i = k; i < n; i++) {
		if (fabs(A[n*i + k]) > m) m = fabs(A[n*i + k]);
	}
	return m;
}

void quad::Lin_eq(int n, double *A, double *b) {
	//direct part
	const double eps = 10e-15;
	const double err = 10e-50;
	double cmax;
	for (int k = 0; k < n - 1; k++) {
		//looking for the non-zero entry
		cmax = col_max(A, n, k);

		if (cmax < err) throw std::invalid_argument("The matrix A is not inversible");
		for (int i = k; i < n; i++) 
			if (fabs(A[n*i + k]) > eps * cmax) {
				swap_rows(A, n, i, k);
				std::swap(b[i], b[k]);
				break;
		}
		//check correctness
		if (fabs(A[n*k + k]) < eps * cmax) throw std::invalid_argument("The matrix A is not inversible");
		//division
		if (fabs(A[n*k + k]) < err) throw std::invalid_argument("The matrix A is not inversible");
		b[k] /= A[n*k+k];
		divide_by(A, n, k, A[n*k + k]);
		//substraction
		for (int i = k + 1; i < n; i++) {
			b[i] -= b[k] * A[n*i + k];
			add_rows(A, n, i, k, -A[n*i + k]);
		}
	}
        //division on the last row
        if (fabs(A[n*(n-1) + n-1]) < err) throw std::invalid_argument("The matrix A is not inversible");
        b[n-1] /= A[n*n - 1];
        divide_by(A, n, n-1, A[n*n - 1]);
	//inverse part
	for (int k = n - 1; k > 0; k--) {
		for (int i = k - 1; i >= 0; i--) {
			b[i] -= b[k] * A[n*i+k];
			add_rows(A, n, i, k, -A[n*i + k]);
		}
	}

}


