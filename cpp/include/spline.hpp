#ifndef SPLINE_HPP
#define SPLINE_HPP

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#include "Eigen/Dense"
#include <vector>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include "omp.h"
#include <chrono>

class StopWatch{
public:
	StopWatch();

	void start();
	void stop();
	void reset();
	void print();
	void print(int cycles);
	double time();

private:
	double time_elapsed;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
};

typedef std::vector<double> Vector;
typedef Eigen::VectorXd EigenVector;
typedef Eigen::ArrayXd EigenArray;

class Spline{
	public:
		// Constructor for cubic spline of f(x)
		Spline(Vector x, Vector f);
		Spline(const EigenArray &x, const EigenArray &f);
		// Copy constructor
		Spline(const Spline& spline);
		// Destructor
		~Spline();

		// Use spline to interpolate f at x = x0
		double evaluate(const double &x0);
		double derivative(const double &x0);
		double derivative2(const double &x0);

		// reconstruct spline on new data set
		void reconstruct(Vector x, Vector f);

	private:
		gsl_spline* _spline;
		gsl_interp_accel* _acc;
};

inline double Spline::evaluate(const double &x0){
	return gsl_spline_eval(_spline, x0, _acc);
}

inline double Spline::derivative(const double &x0){
	return gsl_spline_eval_deriv(_spline, x0, _acc);
}

inline double Spline::derivative2(const double &x0){
	return gsl_spline_eval_deriv2(_spline, x0, _acc);
}

// 2D Interpolation

class Spline2D{
	public:
		// Constructor for cubic spline of f(x)
		Spline2D(Vector x, Vector y, Vector f);
		Spline2D(const EigenArray &x, const EigenArray &y, const EigenArray &f);
		// Copy constructor
		Spline2D(const Spline2D& spline);
		// Destructor
		~Spline2D();

		// Use spline to interpolate f at x = x0
		double evaluate(const double &x0, const double &y0);
		double derivative_x(const double &x0, const double &y0);
		double derivative_y(const double &x0, const double &y0);
		double derivative_xx(const double &x0, const double &y0);
		double derivative_yy(const double &x0, const double &y0);
		double derivative_xy(const double &x0, const double &y0);

		// reconstruct spline on new data set
		void reconstruct(Vector x, Vector y, Vector f);

	private:
		gsl_spline2d* _spline;
		gsl_interp_accel* _xacc;
		gsl_interp_accel* _yacc;
};

inline double Spline2D::evaluate(const double &x0, const double &y0){
	return gsl_spline2d_eval(_spline, y0, x0, _xacc, _yacc);
}

inline double Spline2D::derivative_x(const double &x0, const double &y0){
	return gsl_spline2d_eval_deriv_y(_spline, y0, x0, _xacc, _yacc);
}

inline double Spline2D::derivative_y(const double &x0, const double &y0){
	return gsl_spline2d_eval_deriv_x(_spline, y0, x0, _xacc, _yacc);
}

inline double Spline2D::derivative_xx(const double &x0, const double &y0){
	return gsl_spline2d_eval_deriv_xx(_spline, y0, x0, _xacc, _yacc);
}

inline double Spline2D::derivative_xy(const double &x0, const double &y0){
	return gsl_spline2d_eval_deriv_xy(_spline, y0, x0, _xacc, _yacc);
}

inline double Spline2D::derivative_yy(const double &x0, const double &y0){
	return gsl_spline2d_eval_deriv_yy(_spline, y0, x0, _xacc, _yacc);
}

class Matrix{
public:
	Matrix();
	Matrix(int n);
	Matrix(int n, int m);
	Matrix(int n, int m, Vector A);
	Matrix(int n, int m, double val);

	int rows() const;
	int cols() const;
	int size() const;

	void row_replace(int i, Vector row);
	void col_replace(int j, Vector col);

	Vector row(int i);
	Vector col(int j);

	void reshape(int n, int m);
	Matrix reshaped(int n, int m) const;
	Matrix transpose() const;

	double& operator()(int i, int j);
	const double& operator()(int i, int j) const;

private:
	int _n;
	int _m;
	Vector _A;
};

/////////////////////////////////////////////////////////
////               Basic Interpolators               ////
/////////////////////////////////////////////////////////

class CubicInterpolator{
public:
	CubicInterpolator(double x0, double dx, const Vector &y);
	CubicInterpolator(const Vector &x, const Vector &y);

    CubicInterpolator(double x0, double dx, int nintervals, Matrix cij);
	
    // EigenVector computeDerivatives(double dx, const Vector &y);
    // Eigen::MatrixXd computeDerivativeVector(double dx, const Vector &y);
	double evaluate(const double x);
    double derivative(const double x);
    double derivative2(const double x);

private:
	double evaluateInterval(int i, const double x);
    double evaluateDerivativeInterval(int i, const double x);
    double evaluateSecondDerivativeInterval(int i, const double x);
	void computeSplineCoefficients(double dx, const Vector &y);
	int findInterval(const double x);

	double dx;
	int nintervals;
	double x0;
	Matrix cij;
};

class BicubicInterpolator{
public:
	BicubicInterpolator(const Vector &x, const Vector &y, const Matrix &z);
	BicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Matrix &z);
	double evaluate(const double x, const double y);
    double derivative_x(const double x, const double y);
    double derivative_y(const double x, const double y);
    double derivative_xy(const double x, const double y);
    double derivative_xx(const double x, const double y);
    double derivative_yy(const double x, const double y);
    CubicInterpolator reduce_x(const double x);
    CubicInterpolator reduce_y(const double y);

private:
	double evaluateInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYYInterval(int i, int j, const double x, const double y);
	void computeSplineCoefficients(const Matrix &z);
	int findXInterval(const double x);
	int findYInterval(const double y);

	double dx;
	double dy;
	int nx;
	int ny;
	double x0;
	double y0;
	Eigen::MatrixXd cij;
};

/////////////////////////////////////////////////////////
////            Eigen-based Interpolators            ////
/////////////////////////////////////////////////////////

class EigenCubicInterpolator{
public:
	EigenCubicInterpolator(double x0, double dx, const EigenArray &y);
	EigenCubicInterpolator(const EigenArray &x, const EigenArray &y);

	EigenCubicInterpolator(double x0, double dx, const EigenVector &y);
	EigenCubicInterpolator(const EigenVector &x, const EigenVector &y);

    EigenCubicInterpolator(double x0, double dx, int nintervals, const Eigen::MatrixXd &cij);
	
    EigenVector computeDerivatives(double dx, const EigenVector &y);
    Eigen::MatrixXd computeDerivativeVector(double dx, const EigenVector &y);
	double evaluate(const double x);
    double derivative(const double x);
    double derivative2(const double x);

private:
	double evaluateInterval(int i, const double x);
    double evaluateDerivativeInterval(int i, const double x);
    double evaluateSecondDerivativeInterval(int i, const double x);
	void computeSplineCoefficients(double dx, const EigenArray &y);
	int findInterval(const double x);

	double dx;
	int nintervals;
	double x0;
	Eigen::MatrixXd cij;
};

class EigenBicubicInterpolator{
public:
	EigenBicubicInterpolator(const EigenArray &x, const EigenArray &y, const EigenArray &z);
	EigenBicubicInterpolator(const EigenVector &x, const EigenVector &y, const Eigen::MatrixXd &z);
	EigenBicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Eigen::MatrixXd &z);
    Eigen::MatrixXd computeDerivatives(const Eigen::MatrixXd &z);
    Eigen::MatrixXd computeXDerivatives(const Eigen::MatrixXd &z);
    Eigen::MatrixXd computeYDerivatives(const Eigen::MatrixXd &z);
    Eigen::MatrixXd computeXYDerivatives(const Eigen::MatrixXd &z);
	double evaluate(const double x, const double y);
    double derivative_x(const double x, const double y);
    double derivative_y(const double x, const double y);
    double derivative_xy(const double x, const double y);
    double derivative_xx(const double x, const double y);
    double derivative_yy(const double x, const double y);

	void evaluate(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);
	void derivative_x(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);
	void derivative_y(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);
	void derivative_xy(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);
	void derivative_xx(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);
	void derivative_yy(double z[], const double x[], int nx, const double y[], int ny, int num_threads=0);

	void evaluate(double z[], const double x[], const double y[], int ny, int num_threads=0);
	void derivative_x(double z[], const double x[], const double y[], int ny, int num_threads=0);
	void derivative_y(double z[], const double x[], const double y[], int ny, int num_threads=0);
	void derivative_xy(double z[], const double x[], const double y[], int ny, int num_threads=0);
	void derivative_xx(double z[], const double x[], const double y[], int ny, int num_threads=0);
	void derivative_yy(double z[], const double x[], const double y[], int ny, int num_threads=0);

    EigenCubicInterpolator reduce_x(const double x);
    EigenCubicInterpolator reduce_y(const double y);

private:
	double evaluateInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYYInterval(int i, int j, const double x, const double y);
    Eigen::MatrixXd evaluateInterval(int i, int j, const EigenVector &x, const EigenVector &y);
	void computeSplineCoefficients(const Eigen::MatrixXd &z);
	int findXInterval(const double x);
	int findYInterval(const double y);
    Eigen::ArrayXi findXInterval(const EigenVector &x);
    Eigen::ArrayXi findYInterval(const EigenVector &y);

	double dx;
	double dy;
	int nx;
	int ny;
	double x0;
	double y0;
	Eigen::MatrixXd cij;
};

#endif