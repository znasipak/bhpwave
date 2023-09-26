#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <vector>
#include <algorithm>
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
	void transposeInPlace();

	void set_value(int i, int j, double val);

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

class CubicSpline{
public:
	CubicSpline(double x0, double dx, const Vector &y, int method = 1);
	CubicSpline(const Vector &x, const Vector &y, int method = 1);

    CubicSpline(double x0, double dx, int nintervals, Matrix cij);
	
	double evaluate(const double x);
    double derivative(const double x);
    double derivative2(const double x);

	double getSplineCoefficient(int i, int j);

private:
	double evaluateInterval(int i, const double x);
    double evaluateDerivativeInterval(int i, const double x);
    double evaluateSecondDerivativeInterval(int i, const double x);
	void computeSplineCoefficients(double dx, const Vector &y);
	void computeSplineCoefficientsNaturalFirst(double dx, const Vector &y);
	void computeSplineCoefficientsNotAKnot(double dx, const Vector &y);
	void computeSplineCoefficientsZeroClamped(double dx, const Vector &y);
	void computeSplineCoefficientsE3(double dx, const Vector &y);
	int findInterval(const double x);

	double dx;
	int nintervals;
	double x0;
	Matrix cij;
};

class BicubicSpline{
public:
	BicubicSpline(const Vector &x, const Vector &y, Matrix &z, int method = 3);
	BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, Matrix &z, int method = 3);
	BicubicSpline(const Vector &x, const Vector &y, const Vector &z, int method = 3);
	BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, const Vector &z_vec, int method = 3);
	double evaluate(const double x, const double y);
    double derivative_x(const double x, const double y);
    double derivative_y(const double x, const double y);
    double derivative_xy(const double x, const double y);
    double derivative_xx(const double x, const double y);
    double derivative_yy(const double x, const double y);
    CubicSpline reduce_x(const double x);
    CubicSpline reduce_y(const double y);

private:
	double evaluateInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYYInterval(int i, int j, const double x, const double y);
	Matrix computeSplineCoefficientsDX(Matrix &m_z, int method = 3);
	Matrix computeSplineCoefficientsDY(Matrix &m_z, int method = 3);
	void computeSplineCoefficients(Matrix &z, int method = 3);
	int findXInterval(const double x);
	int findYInterval(const double y);

	double dx;
	double dy;
	int nx;
	int ny;
	double x0;
	double y0;
	Matrix cij;
};

#endif