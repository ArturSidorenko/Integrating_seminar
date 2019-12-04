#pragma once

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<exception>
#include<string>
#include<algorithm>
#include<utility>
#include<cstring>

namespace quad {

	// this thing stores all th coefficients
	//class structure prevents any changes during calculations
	//
	//NB this class stores data for the segment [-1, 1]


	class qu_formula
	{
    private:
		size_t n_;
		double *weights_, *points_;
    public:
		qu_formula();
		qu_formula(size_t n, double *w, double *p) : n_(n), weights_(w), points_(p) {};
        qu_formula(const qu_formula &q); //copy constructor
        qu_formula(qu_formula &&q); //move constructor
        const qu_formula& operator=(const qu_formula &q);
        const qu_formula& operator=(qu_formula &&other);
        ~qu_formula();
        
        int get_n() const {return n_;};
        double get_weight(size_t k) const;
        double get_point(size_t k) const;
	};

	typedef  double(*REAL_FUNC)(double); //the real function

	double simple_quad(REAL_FUNC f, const qu_formula &q, double a, double b); //simple quadrature formula
	double comp_quad(REAL_FUNC f, const qu_formula &q, double a, double b, int N); //compound quadrature formula

	qu_formula Newton_Kotes(); //making quadratures
	qu_formula Gauss();
	qu_formula Naive();

	double col_max(double *A, int n, int k); //to control errors
	void swap_rows(double *A, int n, int p, int q);
	void add_rows(double *A, int n, int p, int q, double m);
	void divide_by(double *A, int n, int p, double d);
	void Lin_eq(int n, double *A, double *b);  //Linear system solving. The matrix A will be spoiled, the answer will be in the vector b

	void find_weights(int n, double *ans, const double *p);
    
}

namespace qu2d {
    //this namespace is for 2d integrating
    typedef std::pair<double, double> point;
    typedef std::pair<size_t, size_t> edge;
    

	typedef double(*REAL_FUNC2D)(double, double); //2d real fucnction

	point get_centre(const point &p1, const point &p2); //finds the centre of the segmet between the given points


	struct side {
		size_t v1, v2, v3; //vertices
		size_t e1, e2, e3; //edges
	};
	struct neib_side {
		int sides;
		size_t first, second;

		neib_side() { sides = 0; }
		void fit(size_t num);
	};




    class partition { //rectangle partition
    private:
        size_t npoints_; //total amount of points
        size_t nedges_; //amount of edges
        size_t nsides_; //amount of sides
        point* points_; //array of points
        edge* edges_; //array of edges
		neib_side* sides_neib_to_edge_; //which sides are neighbouring to the edge
        side* sides_; //array of sides

		bool equal_square_; //whether the square of all the segments is equal
		double sqr_; //square of a triangle 

		double integrate_over_side(REAL_FUNC2D f, const side &k);
    public:
        partition(); //empty constructor
        partition(size_t nx, size_t ny, point down, point upper); //partition of a rectangle
        partition(const std::string &fname); //loader
        ~partition(); //destructor
        
        partition(const partition &other); //copy constructor
        partition(partition &&other); //move onstructor
        
        const partition & operator=(const partition &other);
        const partition & operator=(partition &&other);

		const edge & get_edge(size_t k); //get the edge by its number
		const point & get_point(size_t k); //get the point by its number

		double geron_formula(const side &k);
		double edge_length(const edge &k);

		//integrating
		
		double integrate(REAL_FUNC2D f);


		//saving
		void save(const std::string &fname);
	};
    
}
