#include "Eigen/Dense"
#include <vector>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

typedef std::vector<double> Vector;
typedef Eigen::VectorXd EigenVector;
typedef Eigen::ArrayXd EigenArray;

class Spline{
	public:
		// Constructor for cubic spline of f(x)
		Spline(Vector x, Vector f);
		// Copy constructor
		Spline(const Spline& spline);
		// Destructor
		~Spline();

		// Use spline to interpolate f at x = x0
		double interpolate(const double &x0);
		double derivative(const double &x0);
		double derivative2(const double &x0);

		// reconstruct spline on new data set
		void reconstruct(Vector x, Vector f);

	private:
		gsl_spline* _spline;
		gsl_interp_accel* _acc;
};

inline double Spline::interpolate(const double &x0){
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
		// Copy constructor
		Spline2D(const Spline2D& spline);
		// Destructor
		~Spline2D();

		// Use spline to interpolate f at x = x0
		double interpolate(const double &x0, const double &y0);
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

inline double Spline2D::interpolate(const double &x0, const double &y0){
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

/////////////////////////////////////////////////////////
////            Eigen-based Interpolators            ////
/////////////////////////////////////////////////////////

class CubicInterpolator{
public:
	CubicInterpolator(double x0, double dx, const EigenArray &y);
	CubicInterpolator(const EigenArray &x, const EigenArray &y);

	CubicInterpolator(double x0, double dx, const EigenVector &y);
	CubicInterpolator(const EigenVector &x, const EigenVector &y);

    CubicInterpolator(double dx, const Eigen::MatrixXd &cij);
	
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

class BicubicInterpolator{
public:
	BicubicInterpolator(const EigenArray &x, const EigenArray &y, const EigenArray &z);
	BicubicInterpolator(const EigenVector &x, const EigenVector &y, const Eigen::MatrixXd &z);
	BicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Eigen::MatrixXd &z);
    Eigen::MatrixXd computeDerivatives(Eigen::MatrixXd z);
    Eigen::MatrixXd computeXDerivatives(Eigen::MatrixXd z);
    Eigen::MatrixXd computeYDerivatives(Eigen::MatrixXd z);
    Eigen::MatrixXd computeXYDerivatives(Eigen::MatrixXd z);
	double evaluate(const double x, const double y);
    double derivative_x(const double x, const double y);
    double derivative_y(const double x, const double y);
    double derivative_xy(const double x, const double y);
    double derivative_xx(const double x, const double y);
    double derivative_yy(const double x, const double y);
    Eigen::MatrixXd evaluate(const EigenVector &x, const EigenVector &y);
    CubicInterpolator reduce_x(const double x);
    CubicInterpolator reduce_y(const double y);

private:
	double evaluateInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXYInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeXXInterval(int i, int j, const double x, const double y);
    double evaluateDerivativeYYInterval(int i, int j, const double x, const double y);
    Eigen::MatrixXd evaluateInterval(int i, int j, const EigenVector &x, const EigenVector &y);
	void computeSplineCoefficients(Eigen::MatrixXd z);
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

