#include "integrals.h"

quad::qu_formula::qu_formula() {
    n_ = 0;
    weights_ = nullptr;
    points_ = nullptr;
}

quad::qu_formula::qu_formula(const qu_formula &other) {
    n_ = other.n_;
    weights_ = new double[n_];
    points_ = new double[n_];
    std::memcpy(weights_, other.weights_, n_);
    std::memcpy(points_, other.points_, n_);
}

quad::qu_formula::~qu_formula() {
    if(points_  != nullptr) delete [] points_;
	if(weights_ != nullptr) delete [] weights_;
}

const quad::qu_formula& quad::qu_formula::operator=(const qu_formula &other) {
    if(this != &other) {
        delete [] points_;
        delete [] weights_;
        n_ = other.n_;
        weights_ = new double[n_];
        points_ = new double[n_];
        std::memcpy(weights_, other.weights_, n_*sizeof(double));
        std::memcpy(points_, other.points_, n_ * sizeof(double));
    }
    return *this;
}

quad::qu_formula::qu_formula(qu_formula &&other) {
    weights_ = other.weights_;
    points_ = other.points_;
    std::swap(n_, other.n_);
    other.weights_ = nullptr;
    other.points_ = nullptr;

}

const quad::qu_formula& quad::qu_formula::operator=(qu_formula &&other) {
    if(this !=&other) {
        delete [] points_;
        delete [] weights_;
        n_ = other.n_;
        weights_ = other.weights_;
        points_ = other.points_;
        other.weights_ = nullptr;
        other.points_ = nullptr;
    }
    return *this;
}

double quad::qu_formula::get_point(size_t k) const{
    if((k>=n_)) throw std::out_of_range("The required point index is out of range");
    return points_[k];
}

double quad::qu_formula::get_weight(size_t k) const{
    if((k>=n_)) throw std::out_of_range("The required weight index is out of range");
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

		if (cmax < err) throw std::invalid_argument("The matrix A is not invertible");
		for (int i = k; i < n; i++) 
			if (fabs(A[n*i + k]) > eps * cmax) {
				swap_rows(A, n, i, k);
				std::swap(b[i], b[k]);
				break;
		}
		//check correctness
		if (fabs(A[n*k + k]) < eps * cmax) throw std::invalid_argument("The matrix A is not invertible");
		//division
		if (fabs(A[n*k + k]) < err) throw std::invalid_argument("The matrix A is not invertible");
		b[k] /= A[n*k+k];
		divide_by(A, n, k, A[n*k + k]);
		//substraction
		for (int i = k + 1; i < n; i++) {
			b[i] -= b[k] * A[n*i + k];
			add_rows(A, n, i, k, -A[n*i + k]);
		}
	}
        //division on the last row
        if (fabs(A[n*(n-1) + n-1]) < err) throw std::invalid_argument("The matrix A is not invertible");
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


qu2d::partition::partition() {
	npoints_ = nedges_ = nsides_ = 0;
	points_ = nullptr;
	edges_ = nullptr;
	sides_ = nullptr;
	sides_neib_to_edge_ = nullptr;
}

qu2d::partition::~partition() {
	if (points_ != nullptr) delete[] points_;
	if (edges_ != nullptr) delete[] edges_;
	if (sides_ != nullptr) delete[] sides_;
	if (sides_neib_to_edge_ != nullptr) delete[] sides_neib_to_edge_;
}

qu2d::partition::partition(const partition &other): 
	npoints_(other.npoints_),
	nedges_(other.nedges_), nsides_(other.nsides_), sqr_(other.sqr_), equal_square_(other.equal_square_)
{
	points_ = new point[npoints_];
	edges_ = new edge[nedges_];
	sides_ = new side[nsides_];
	sides_neib_to_edge_ = new neib_side[nedges_];

	memcpy(points_, other.points_, npoints_ * sizeof(point));
	memcpy(edges_, other.edges_, nedges_ * sizeof(edge));
	memcpy(sides_, other.sides_, nsides_ * sizeof(side));
	memcpy(sides_neib_to_edge_, other.sides_neib_to_edge_, nedges_ * sizeof(neib_side));

}

qu2d::partition::partition(partition &&other) :
	npoints_(other.npoints_),
	nedges_(other.nedges_), nsides_(other.nsides_), sqr_(other.sqr_), equal_square_(other.equal_square_)
{
	points_ = other.points_;
	edges_ = other.edges_;
	sides_ = other.sides_;
	sides_neib_to_edge_ = other.sides_neib_to_edge_;

	other.points_ = nullptr;
	other.edges_ = nullptr;
	other.sides_ = nullptr;
	other.sides_neib_to_edge_ = nullptr;

}

const qu2d::partition & qu2d::partition::operator=(const qu2d::partition &other) {
	if (this != &other) {
		npoints_ = other.npoints_;
		nedges_ = other.nedges_;
		nsides_ = other.nsides_;
		sqr_ = other.sqr_;
		equal_square_ = other.equal_square_;

		delete[] points_;
		delete[] edges_;
		delete[] sides_;
		delete[] sides_neib_to_edge_;

		points_ = new point[npoints_];
		edges_ = new edge[nedges_];
		sides_ = new side[nsides_];
		sides_neib_to_edge_ = new neib_side[nedges_];

		memcpy(points_, other.points_, npoints_ * sizeof(point));
		memcpy(edges_, other.edges_, nedges_ * sizeof(edge));
		memcpy(sides_, other.sides_, nsides_ * sizeof(side));
		memcpy(sides_neib_to_edge_, other.sides_neib_to_edge_, nedges_ * sizeof(neib_side));
	}
	return *this;
}

const qu2d::partition & qu2d::partition::operator=(qu2d::partition &&other) {
	if (this != &other) {
		npoints_ = other.npoints_;
		nedges_ = other.nedges_;
		nsides_ = other.nsides_;
		sqr_ = other.sqr_;
		equal_square_ = other.equal_square_;

		points_ = other.points_;
		edges_ = other.edges_;
		sides_ = other.sides_;
		sides_neib_to_edge_ = other.sides_neib_to_edge_;

		other.points_ = nullptr;
		other.edges_ = nullptr;
		other.sides_ = nullptr;
		other.sides_neib_to_edge_ = nullptr;

	}

	return *this;
}

//this is the most interesting thing 
qu2d::partition::partition(size_t nx, size_t ny, point down, point upper)
{
	//calculate numbers
	npoints_ = (nx + 1)*(ny + 1);
	size_t vert_edges = (nx + 1)*ny;
	size_t horis_edges = (ny + 1)*nx;
	size_t diag_edges = nx * ny;
	nedges_ = horis_edges + vert_edges + diag_edges;
	nsides_ = nx * ny * 2; 

	//assign memory
	points_ = new point[npoints_];
	edges_ = new edge[nedges_];
	sides_ = new side[nsides_];
	sides_neib_to_edge_ = new neib_side[nedges_];

	//find out vertices
	double dx = (upper.first - down.first) / nx;
	double dy = (upper.second - down.second) / ny;

	sqr_ = 0.5 * dx * dy;
	equal_square_ = true;

	for (size_t i = 0; i <= nx; i++) {
		for (size_t j = 0; j <= ny; j++) {
			points_[(ny+1)*i + j] = point(down.first + i * dx, down.second + j * dy);
		}
	}

	//making vertical edges 
	for (size_t i = 0; i <= nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			edges_[ny*i + j] = edge((ny+1)*i + j, (ny + 1)*i + j+1);
		}
	}
	//making horizontal edges
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j <= ny; j++) {
			edges_[vert_edges + (ny+1)*i + j] = edge((ny + 1)*i + j, (ny + 1)*(i + 1) + j);
		}
	}
	//making diagonal edges (rising)
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			edges_[vert_edges + horis_edges + ny *i + j] = edge((ny + 1)*i + j, (ny + 1)*(i + 1) + j + 1);
		}
	}


	//making sides
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			//the down triangle
			sides_[2 * (ny*i + j)].v1 = (ny + 1) * i + j;
			sides_[2 * (ny*i + j)].v2 = (ny + 1) * i + j + 1;
			sides_[2 * (ny*i + j)].v3 = (ny + 1) * (i + 1) + j;
			sides_[2 * (ny*i + j)].e1 = ny * i + j; //vertical edge
			sides_[2 * (ny*i + j)].e2 = vert_edges + (ny + 1)*i + j; //gorisontal edge
			sides_[2 * (ny*i + j)].e3 = vert_edges + horis_edges + ny * i + j; //diagonal edge
			//neighbouring sides
			sides_neib_to_edge_[ny * i + j].fit(2 * (ny*i + j));
			sides_neib_to_edge_[vert_edges + (ny + 1)*i + j].fit(2 * (ny*i + j));
			sides_neib_to_edge_[vert_edges + horis_edges + ny * i + j].fit(2 * (ny*i + j));

			//the upper triangle
			sides_[2 * (ny*i + j) + 1].v1 = (ny + 1) * i + 1 + j + 1;
			sides_[2 * (ny*i + j) + 1].v2 = (ny + 1) * i + j + 1;
			sides_[2 * (ny*i + j) + 1].v3 = (ny + 1) * (i + 1) + j;
			sides_[2 * (ny*i + j) + 1].e1 = ny * (i + 1) + j; //vertical edge
			sides_[2 * (ny*i + j) + 1].e2 = vert_edges + (ny + 1)*i + j + 1; //gorisontal edge
			sides_[2 * (ny*i + j) + 1].e3 = vert_edges + horis_edges + ny * i + j; //diagonal edge

			//neibouring sides
			sides_neib_to_edge_[ny * (i + 1) + j].fit(2 * (ny*i + j) + 1);
			sides_neib_to_edge_[vert_edges + (ny + 1)*i + j + 1].fit(2 * (ny*i + j) + 1);
			sides_neib_to_edge_[vert_edges + horis_edges + ny * i + j].fit(2 * (ny*i + j) + 1);

		}
	}

}

void qu2d::neib_side::fit(size_t num) {
	switch (sides)
	{
	case 0:
		first = num;
		break;
	case 1:
		second = num;
		break;
	default:
		throw(std::out_of_range("Neib_size filling error occured\n"));
		break;
	}
	sides++;
}

const qu2d::edge & qu2d::partition::get_edge(size_t k) {
	if (k >= nedges_) throw std::out_of_range("The required edge index is out of range ");
	return edges_[k];
}

const qu2d::point & qu2d::partition::get_point(size_t k) {
	if (k >= npoints_) throw std::out_of_range("The required point index is out of range ");
	return points_[k];
}

qu2d::point qu2d::get_centre(const point &p1, const point &p2) {
	double x, y;
	x = 0.5*(p1.first + p2.first);
	y = 0.5*(p1.second + p2.second);
	return point(x, y);
}

double qu2d::partition::integrate_over_side(qu2d::REAL_FUNC2D f, const side &k) {
	const point& p1 = points_[k.v1];
	const point& p2 = points_[k.v2];
	const point& p3 = points_[k.v3];

	point m1 = get_centre(p1, p2);
	point m2 = get_centre(p2, p3);
	point m3 = get_centre(p1, p3);

	double z1 = f(p1.first, p1.second);
	double z2 = f(p2.first, p2.second);
	double z3 = f(p3.first, p3.second);

	double sqr = equal_square_ ? sqr_ : geron_formula(k);

	return sqr * (z1 + z2 + z3) / 3;

}

double qu2d::partition::integrate(qu2d::REAL_FUNC2D f) {
	double ans = 0;
	for (size_t i = 0; i < nsides_; i++) {
		ans += integrate_over_side(f, sides_[i]);
	}
	return ans;
}

void qu2d::partition::save(const std::string &fname) {
	std::ofstream out(fname);

	//writing key numbers
	out << npoints_ << " "<< nedges_ << " " << nsides_ << "\n";
	out << (int)(equal_square_) << "\n";

	out << "\n";

	//writing point coordinates
	for (size_t i = 0; i < npoints_; i++) {
		out << points_[i].first << " " << points_[i].second << "\n";
	}

	out << "\n";

	//writing edges
	for (size_t i = 0; i < nedges_; i++) {
		out << edges_[i].first << " " << edges_[i].second<<"\n";
	}
	out << "\n";
	//writing sides
	for (size_t i = 0; i < nsides_; i++) {
		out << sides_[i].v1 << " " << sides_[i].v2 << " " << sides_[i].v3 << " " << 
			sides_[i].e1 << " " << sides_[i].e2 << " " << sides_[i].e3 << "\n";
	}
	out << "\n";
	//writing edge neibours
	for (size_t i = 0; i < nedges_; i++) {
		out << sides_neib_to_edge_[i].sides << " ";
		if (sides_neib_to_edge_[i].sides >= 1) out << sides_neib_to_edge_[i].first << " ";
		if (sides_neib_to_edge_[i].sides >= 2) out << sides_neib_to_edge_[i].second;
		out << "\n";
	}
	out.close();
}

qu2d::partition::partition(const std::string &fname) {
	std::ifstream in(fname);

	//key numbers
	in >> npoints_;
	in >> nedges_;
	in >> nsides_;
	in >> equal_square_;

	//memory assignment
	points_ = new point[npoints_];
	edges_ = new edge[nedges_];
	sides_ = new side[nsides_];
	sides_neib_to_edge_ = new neib_side[nedges_];
	
	//point coordinates
	for (size_t i = 0; i < npoints_; i++) {
		in >> points_[i].first;
		in >> points_[i].second;
	}

	//edge coordinates
	for (size_t i = 0; i < nedges_; i++) {
		in >> edges_[i].first;
		in >> edges_[i].second;
	}

	//sides
	for (size_t i = 0; i < nsides_; i++) {
		in >> sides_[i].v1;
		in >> sides_[i].v2;
		in >> sides_[i].v3;
		in >> sides_[i].e1;
		in >> sides_[i].e2;
		in >> sides_[i].e3;
	}

	//side neighbours
	for (size_t i = 0; i < nedges_; i++) {
		in >> sides_neib_to_edge_[i].sides;
		if (sides_neib_to_edge_[i].sides >= 1) in >> sides_neib_to_edge_[i].first;
		if (sides_neib_to_edge_[i].sides >= 2) in >> sides_neib_to_edge_[i].second;
	}

	if (equal_square_) sqr_ = geron_formula(sides_[0]);
	in.close();
}

double qu2d::partition::geron_formula(const side &k) {
	double a = edge_length(edges_[k.e1]), b = edge_length(edges_[k.e2]), c = edge_length(edges_[k.e3]);

	double p = 0.5*(a + b + c);

	return sqrt(p * (p - a) * (p - b) * (p - c));
}

double qu2d::partition::edge_length(const edge &k) {
	double x = points_[k.second].first - points_[k.first].first;
	double y = points_[k.second].second - points_[k.first].second;

	return sqrt(x * x + y * y);
	
}