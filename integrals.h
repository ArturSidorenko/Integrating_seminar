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
    typedef std::pair<size_t, size_t> neib_side;
    typedef std::triple<size_t, size_t, size_t> side;
    class partition { //rectangle partition
    private:
        size_t nx_, ny_; //amount of pieces along each axis
        size_t npoints_; //total amount of points
        size_t nedge_; //amount of edges
        size_t nsides_; //amount of sides
        point* points_; //array of points
        edge* edges_; //array of edges
        side* sides_; //array of sides
        
    public:
        partition(); //empty constructor
        partition(size_t nx, size_t ny, point down, point upper); //partition of a rectangle
        partition(std::string &fname); //loader
        ~partition(); //destructor
        
        partition(const partition &other); //copy constructor
        partition(partition &&other); //move onstructor
        
        const parition & operator=(const partition &other);
        const parition & operator=(partition &&other);
    }
    
    class partition_iterator: public std::iterator<std::input_iterator_tag, partition> {
        friend class partition;
    private:
        partition_iterator(partition *p);
    public:
        partition_iterator(const partition_iterator &it);
    }
    
}
