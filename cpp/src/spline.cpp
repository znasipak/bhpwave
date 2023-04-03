#include "spline.hpp"
#include <iostream>

#define ENDPOINT_TOL 1.e-10

Spline::Spline(Vector x, Vector f){
	_acc = gsl_interp_accel_alloc();
	_spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
	if(x[0] > x[1]){
		std::reverse(x.begin(), x.end());
		std::reverse(f.begin(), f.end());
	}
	gsl_spline_init(_spline, &x[0], &f[0], x.size());
}

Spline::Spline(const Spline& obj){
	_spline = new gsl_spline;
	_acc = new gsl_interp_accel;

	*_spline = *obj._spline;
	*_acc = *obj._acc;
}

Spline::~Spline(){
	delete(_acc);
	delete(_spline);
}

void Spline::reconstruct(Vector x, Vector f){
	delete(_acc);
	delete(_spline);

	_acc = gsl_interp_accel_alloc();
	_spline = gsl_spline_alloc(gsl_interp_cspline, x.size());

	if(x[0] > x[1]){
		std::reverse(x.begin(), x.end());
		std::reverse(f.begin(), f.end());
	}
	gsl_spline_init(_spline, &x[0], &f[0], x.size());
}

// 2D Interpolation

// x and y are switched because GSL's ordering for f is the transpose of how I typically order my arrays
// in other words,
// GSL ordering is f_ij = f[j*xsize + i]
// my ordering is f_ij = f[i*ysize + j]
// for i = 0, ..., xsize; j = 0, ..., ysize
Spline2D::Spline2D(Vector y, Vector x, Vector f){
	if(x.size()*y.size() != f.size()){
		std::cout << "(SPLINE) Error: (x, y) array with dimensions (" << x.size() << ", " << y.size() << ") not compabitble with f(x, y) array of length " << f.size() << "\n";
	}

	_xacc = gsl_interp_accel_alloc();
	_yacc = gsl_interp_accel_alloc();
	_spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, x.size(), y.size());
	Vector fnew(f.size());
	if(x[0] > x[1] && y[0] > y[1]){
		std::reverse(x.begin(), x.end());
		std::reverse(y.begin(), y.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + x.size() - i - 1];
			}
		}
	}else if(x[0] > x[1]){
		std::reverse(x.begin(), x.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[j*x.size() + x.size() - i - 1];
			}
		}
	}else if(y[0] > y[1]){
		std::reverse(y.begin(), y.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + i];
			}
		}
	}else{
		for(size_t j = 0; j < f.size(); j++){
			fnew[j] = f[j];
		}
	}

	// for(size_t i = 1; i < x.size(); i++){
	// 	if(x[i] <= x[i-1]){
	// 		std::cout << "x["<<i<<"] = "<<x[i]<<"\n";
	// 	}
	// }
	//
	// std::cout << "Foo\n";
	// for(size_t i = 1; i < y.size(); i++){
	// 	if(y[i] <= y[i-1]){
	// 		std::cout << "y["<<i<<"] = "<<y[i]<<"\n";
	// 	}
	// }
	gsl_spline2d_init(_spline, &x[0], &y[0], &fnew[0], x.size(), y.size());
}

// Spline2D::Spline2D(const EigenArray &yarray, const EigenArray &xarray, const EigenArray &farray){
// 	Vector x(xarray.size());
// 	Vector y(yarray.size());
// 	Vector f(farray.size());
// 	Eigen::ArrayXd::Map(&x[0], xarray.size()) = xarray;
// 	Eigen::ArrayXd::Map(&y[0], yarray.size()) = yarray;
// 	Eigen::ArrayXd::Map(&f[0], farray.size()) = farray;

// 	if(x.size()*y.size() != f.size()){
// 		std::cout << "(SPLINE) Error: (x, y) array with dimensions (" << x.size() << ", " << y.size() << ") not compabitble with f(x, y) array of length " << f.size() << "\n";
// 	}

// 	_xacc = gsl_interp_accel_alloc();
// 	_yacc = gsl_interp_accel_alloc();
// 	_spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, x.size(), y.size());
// 	Vector fnew(f.size());
// 	if(x[0] > x[1] && y[0] > y[1]){
// 		std::reverse(x.begin(), x.end());
// 		std::reverse(y.begin(), y.end());
// 		for(size_t i = 0; i < x.size(); i++){
// 			for(size_t j = 0; j < y.size(); j++){
// 				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + x.size() - i - 1];
// 			}
// 		}
// 	}else if(x[0] > x[1]){
// 		std::reverse(x.begin(), x.end());
// 		for(size_t i = 0; i < x.size(); i++){
// 			for(size_t j = 0; j < y.size(); j++){
// 				fnew[j*x.size() + i] = f[j*x.size() + x.size() - i - 1];
// 			}
// 		}
// 	}else if(y[0] > y[1]){
// 		std::reverse(y.begin(), y.end());
// 		for(size_t i = 0; i < x.size(); i++){
// 			for(size_t j = 0; j < y.size(); j++){
// 				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + i];
// 			}
// 		}
// 	}else{
// 		for(size_t j = 0; j < f.size(); j++){
// 			fnew[j] = f[j];
// 		}
// 	}
// 	gsl_spline2d_init(_spline, &x[0], &y[0], &fnew[0], x.size(), y.size());
// }

Spline2D::Spline2D(const Spline2D& obj){
	_spline = new gsl_spline2d;
	_xacc = new gsl_interp_accel;
	_yacc = new gsl_interp_accel;

	*_spline = *obj._spline;
	*_xacc = *obj._xacc;
	*_yacc = *obj._yacc;
}

Spline2D::~Spline2D(){
	delete(_xacc);
	delete(_yacc);
	delete(_spline);
}

// see above for ordering of arguments
void Spline2D::reconstruct(Vector y, Vector x, Vector f){
	delete(_xacc);
	delete(_yacc);
	delete(_spline);

	_xacc = gsl_interp_accel_alloc();
	_yacc = gsl_interp_accel_alloc();
	_spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, x.size(), y.size());

	Vector fnew(f.size());
	if(x[0] > x[1] && y[0] > y[1]){
		std::reverse(x.begin(), x.end());
		std::reverse(y.begin(), y.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + x.size() - i - 1];
			}
		}
	}else if(x[0] > x[1]){
		std::reverse(x.begin(), x.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[(y.size() - j - 1)*x.size() + i];
			}
		}
	}else if(y[0] > y[1]){
		std::reverse(y.begin(), y.end());
		for(size_t i = 0; i < x.size(); i++){
			for(size_t j = 0; j < y.size(); j++){
				fnew[j*x.size() + i] = f[j*x.size() + y.size() - i - 1];
			}
		}
	}else{
		for(size_t j = 0; j < f.size(); j++){
			fnew[j] = f[j];
		}
	}
	gsl_spline2d_init(_spline, &x[0], &y[0], &fnew[0], x.size(), y.size());
}

///////////////////////////////////////////////////////////////////
//////////////              Matrix Class           ////////////////
///////////////////////////////////////////////////////////////////

Matrix::Matrix() {}
Matrix::Matrix(int n): _n(n), _m(n), _A(n*n) {}
Matrix::Matrix(int n, int m): _n(n), _m(m), _A(n*m) {}
Matrix::Matrix(int n, int m, Vector A): _n(n), _m(m), _A(n*m) {
	if(A.size() == _A.size()){
		_A = A;
	}
}
Matrix::Matrix(int n, int m, double val): _n(n), _m(m), _A(n*m, val) {}

int Matrix::rows() const{
	return _n;
}
int Matrix::cols() const{
	return _m;
}
int Matrix::size() const{
	return _A.size();
}

void Matrix::row_replace(int i, Vector row){
	for(int j = 0; j < _m; j++){
		_A[j*_m + i] = row[j];
	}
}
void Matrix::col_replace(int i, Vector col){
	for(int j = 0; j < _n; j++){
		_A[i*_m + j] = col[j];
	}
}

Vector Matrix::row(int i){
	Vector row(_m);
	for(int j = 0; j < _m; j++){
		row[j] = _A[j*_m + i];
	}
	return row;
}
Vector Matrix::col(int i){
	Vector col(_n);
	for(int j = 0; j < _n; j++){
		col[j] = _A[i*_m + j];
	}
	return col;
}

void Matrix::reshape(int n, int m){
	_n = n;
	_m = m;
}

Matrix Matrix::reshaped(int n, int m) const{
	return Matrix(n, m , _A);
}

Matrix Matrix::transpose() const{
	Matrix AT(_m, _n);
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < _n; i++){
				for(int j = 0; j < _m; j++){
					AT(j, i) = _A[i*_m + j];
				}
			}
	}
	return AT;
}

void Matrix::transposeInPlace(){
	Vector AT(_A.size());
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for(int i = 0; i < _n; i++){
			for(int j = 0; j < _m; j++){
				AT[j*_n + i] = _A[i*_m + j];
			}
		}

		#pragma omp for
		for(int i = 0; i < _n*_m; i++){
			_A[i] = AT[i];
		}
	}
	int m = _n;
	_n = _m;
	_m = m;
}

double& Matrix::operator()(int i, int j){
	return _A[i*_m + j];
}
const double& Matrix::operator()(int i, int j) const{
	return _A[i*_m + j];
}

StopWatch::StopWatch():time_elapsed(0.), t1(std::chrono::high_resolution_clock::now()), t2(t1) {}
void StopWatch::start(){
	t1 = std::chrono::high_resolution_clock::now();
	t2 = t1;
}

void StopWatch::stop(){
	t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	time_elapsed += time_span.count();
}

void StopWatch::reset(){
	time_elapsed = 0.;
}

void StopWatch::print(){
	std::cout << "It took me "<< time_elapsed <<" seconds total.";
	std::cout << std::endl;
}

void StopWatch::print(int cycles){
	std::cout << "It took me "<< time_elapsed/cycles <<" seconds per cycle, "<< time_elapsed <<" seconds total.";
	std::cout << std::endl;
}

double StopWatch::time(){
	return time_elapsed;
}

//////////////////////////////////////////////////////////////////
//////////////          BicubicInterpolator       ////////////////
//////////////////////////////////////////////////////////////////

BicubicInterpolator::BicubicInterpolator(const Vector &x, const Vector &y, const Matrix &z): BicubicInterpolator(x[0], x[1] - x[0], x.size() - 1, y[0], y[1] - y[0], y.size() - 1, z) {}
BicubicInterpolator::BicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Matrix &z): dx(dx), dy(dy), nx(nx), ny(ny), x0(x0), y0(y0), cij(nx, 16*ny) {
	if(nx + 1 != z.rows() && ny + 1 != z.cols()){
		if(nx + 1 == z.cols() && ny + 1 == z.rows()){
			// switch x and y
			cij.transposeInPlace();
			computeSplineCoefficients(z);
		}else if((nx + 1)*(ny + 1) == z.size()){
			Matrix m_z = z.reshaped(ny + 1, nx + 1).transpose();
			computeSplineCoefficients(m_z);
		}else{
			std::cout << "ERROR: Indices of vectors and matrices do not match \n";
		}
	}else{
		computeSplineCoefficients(z);
	}
}

double BicubicInterpolator::evaluate(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateInterval(i, j, x, y);
}

double BicubicInterpolator::derivative_x(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXInterval(i, j, x, y)/dx;
}

double BicubicInterpolator::derivative_y(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeYInterval(i, j, x, y)/dy;
}

double BicubicInterpolator::derivative_xy(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXYInterval(i, j, x, y)/dx/dy;
}

double BicubicInterpolator::derivative_xx(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXXInterval(i, j, x, y)/dx/dx;
}

double BicubicInterpolator::derivative_yy(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeYYInterval(i, j, x, y)/dy/dy;
}


void BicubicInterpolator::computeSplineCoefficients(const Matrix &m_z){
	StopWatch watch;

	watch.start();
	Matrix lmat(4, 4, 0.);
	lmat(0, 0) = 1.;
	lmat(1, 2) = 1.;
	lmat(2, 0) = -3.;
	lmat(2, 1) = 3.;
	lmat(2, 2) = -2.;
	lmat(2, 3) = -1.;
	lmat(3, 0) = 2.;
	lmat(3, 1) = -2.;
	lmat(3, 2) = 1.;
	lmat(3, 3) = 1.;
	
	// 4th-order central difference coefficients
	Vector diffopCentral = {1./12., -8./12., 0., 8./12., -1./12.};
	// 4th-order forward difference coefficients
	Vector diffopForward = {-25./12., 4., -3., 4./3., -1./4};
	Matrix m_zdy(m_z.rows(), m_z.cols());

	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdy.rows(); i++){
				for(int j = 0; j < 2; j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_z(i, j + k)*diffopForward[k];
					}
					m_zdy(i, j) = sum;
				}
			}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdy.rows(); i++){
				for(int j = 2; j < m_zdy.cols() - 2; j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_z(i, j - 2 + k)*diffopCentral[k];
					}
					m_zdy(i, j) = sum;
				}
			}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdy.rows(); i++){
				for(int j = m_zdy.cols() - 2; j < m_zdy.cols(); j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += -m_z(i, j - k)*diffopForward[k];
					}
					m_zdy(i, j) = sum;
				}
			}
	}
	watch.stop();
	watch.print();
	watch.reset();

	watch.start();
	// Might be more efficient to transpose m_z here to take advantage of
	// row-major ordering
	Matrix m_zT = m_z.transpose();
	Matrix m_zdxT(m_zT.rows(), m_zT.cols());

	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for(int i = 0; i < m_zdxT.rows(); i++){
			for(int j = 0; j < 2; j++){
				double sum = 0.;
				for(int k = 0; k < 5; k++){
					sum += m_zT(i, j + k)*diffopForward[k];
				}
				m_zdxT(i, j) = sum;
			}
		}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdxT.rows(); i++){
				for(int j = 2; j < m_zdxT.cols() - 2; j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_zT(i, j - 2 + k)*diffopCentral[k];
					}
					m_zdxT(i, j) = sum;
				}
			}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdxT.rows(); i++){
				for(int j = m_zdxT.cols() - 2; j < m_zdxT.cols(); j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += -m_zT(i, j - k)*diffopForward[k];
					}
					m_zdxT(i, j) = sum;
				}
			}
	}
	Matrix m_zdx = m_zdxT.transpose();
	watch.stop();
	watch.print();
	watch.reset();
	
	watch.start();
	Matrix m_zdxdy(m_zdx.rows(), m_zdx.cols());
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdxdy.rows(); i++){
				for(int j = 0; j < 2; j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_zdx(i, j + k)*diffopForward[k];
					}
					m_zdxdy(i, j) = sum;
				}
			}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdxdy.rows(); i++){
				for(int j = 2; j < m_zdxdy.cols() - 2; j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_zdx(i, j - 2 + k)*diffopCentral[k];
					}
					m_zdxdy(i, j) = sum;
				}
			}
	}
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
			for(int i = 0; i < m_zdxdy.rows(); i++){
				for(int j = m_zdxdy.cols() - 2; j < m_zdxdy.cols(); j++){
					double sum = 0.;
					for(int k = 0; k < 5; k++){
						sum += m_zdx(i, j - k)*diffopForward[k];
					}
					m_zdxdy(i, j) = sum;
				}
			}
	}
	watch.stop();
	watch.print();
	watch.reset();

	// now this part we just have to accept as being inefficient because we
	// are mixing rows and columns no matter what. The important thing is that
	// we will store the relevant cofficients close to one another in memory

	watch.start();
	#pragma omp parallel
	{
		#pragma omp for collapse(2) schedule(dynamic, 16)
			for(int i = 0; i < nx; i++){
				for(int j = 0; j < ny; j++){
					Matrix dmat(4, 4);
					dmat(0, 0) = m_z(i, j); // f(0,0)
					dmat(0, 1) = m_z(i, j + 1); // f(0,1)
					dmat(0, 2) = m_zdy(i, j); // fy(0,0)
					dmat(0, 3) = m_zdy(i, j + 1); // fy(0,1)
					dmat(1, 0) = m_z(i + 1, j); // f(1,0)
					dmat(1, 1) = m_z(i + 1, j + 1); // f(1,1)
					dmat(1, 2) = m_zdy(i + 1, j); // fy(1,0)
					dmat(1, 3) = m_zdy(i + 1, j + 1); // fy(1,1)
					dmat(2, 0) = m_zdx(i, j); // fx(0,0)
					dmat(2, 1) = m_zdx(i, j + 1); // fx(0,1)
					dmat(2, 2) = m_zdxdy(i, j); // fxy(0,0)
					dmat(2, 3) = m_zdxdy(i, j + 1); // fxy(0,1)
					dmat(3, 0) = m_zdx(i + 1, j); // fx(1,0)
					dmat(3, 1) = m_zdx(i + 1, j + 1); // fx(1,1)
					dmat(3, 2) = m_zdxdy(i + 1, j); // fxy(1,0)
					dmat(3, 3) = m_zdxdy(i + 1, j + 1); // fxy(1,1)

					// this part is slow. Just lots of matrix multiplication
					Matrix Dmat(4, 4);
					for(int k = 0; k < 4; k++){
						for(int l = 0; l < 4; l++){
							for(int m = 0; m < 4; m++){
								Dmat(l, k) += dmat(k, m)*lmat(k, m);
							}
							for(int m = 0; m < 4; m++){
								cij(i, 16*j + 4*k + l) += lmat(k, m)*Dmat(k, m);
							}
						}
					}
					// for(int k = 0; k < 4; k++){
					// 	for(int l = 0; l < 4; l++){
					// 		for(int m = 0; m < 4; m++){
					// 			cij(i, 16*j + 4*k + l) += lmat(k, m)*Dmat(k, m);
					// 		}
					// 	}
					// }
				}
			}
	}
	watch.stop();
	watch.print();
	watch.reset();
}

double BicubicInterpolator::evaluateInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	double zvec[4] = {0., 0., 0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 4; k++){
		for(int l = 0; l < 4; l++){
			zvec[k] += cij(i, 16*j + 4*k + l)*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 4; k++){
		result += xvec[k]*zvec[k];
	}

	return result;
}

double BicubicInterpolator::evaluateDerivativeXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	double xvec[3] = {1., 2.*xbar, 3.*xbar*xbar};
	double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	double zvec[3] = {0., 0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 3; k++){
		for(int l = 0; l < 4; l++){
			zvec[k] += cij(i, 16*j + 4*(1 + k) + l)*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 3; k++){
		result += xvec[k]*zvec[k];
	}

	return result;
}

double BicubicInterpolator::evaluateDerivativeYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;

	double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	double yvec[3] = {1., 2.*ybar, 3.*ybar*ybar};
	double zvec[4] = {0., 0., 0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 4; k++){
		for(int l = 0; l < 3; l++){
			zvec[k] += cij(i, 16*j + 4*k + (l + 1))*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 4; k++){
		result += xvec[k]*zvec[k];
	}
	
	return result;
}

double BicubicInterpolator::evaluateDerivativeXYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	double xvec[3] = {1., 2.*xbar, 3.*xbar*xbar};
	double yvec[3] = {1., 2.*ybar, 3.*ybar*ybar};
	double zvec[3] = {0., 0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 3; k++){
		for(int l = 0; l < 3; l++){
			zvec[k] += cij(i, 16*j + 4*(k + 1) + (l + 1))*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 3; k++){
		result += xvec[k]*zvec[k];
	}
	
	return result;
}

double BicubicInterpolator::evaluateDerivativeXXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	double xvec[2] = {2., 6.*xbar};
	double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	double zvec[2] = {0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 2; k++){
		for(int l = 0; l < 4; l++){
			zvec[k] += cij(i, 16*j + 4*(2 + k) + l)*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 2; k++){
		result += xvec[k]*zvec[k];
	}

	return result;
}

double BicubicInterpolator::evaluateDerivativeYYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	double yvec[2] = {2., 6.*ybar};
	double zvec[4] = {0., 0., 0., 0.};
	double result = 0.;

	// zj = cij*yi
	for(int k = 0; k < 4; k++){
		for(int l = 0; l < 2; l++){
			zvec[k] += cij(i, 16*j + 4*k + (l + 2))*yvec[l];
		}
	}

	// result = xj*zj
	for(int k = 0; k < 4; k++){
		result += xvec[k]*zvec[k];
	}
	
	return result;
}

int BicubicInterpolator::findXInterval(const double x){
	int i = static_cast<int>((x-x0)/dx);
    if(i >= nx){
        return nx - 1;
    }
	if( i < 0){
		return 0;
	}
	return i;
}

int BicubicInterpolator::findYInterval(const double y){
	int i = static_cast<int>((y-y0)/dy);
    if(i >= ny){
        return ny - 1;
    }
	if( i < 0){
		return 0;
	}
	return i;
}

CubicInterpolator BicubicInterpolator::reduce_x(const double x){
    int i = findXInterval(x);
    double xbar = (x - x0 - i*dx)/dx;

	Matrix cubicCij(ny, 4);
	double xvec[4] = {1., xbar, xbar*xbar, xbar*xbar*xbar};

	// zj = xi*cij
	for(int j = 0; j < ny; j++){
		for(int k = 0; k < 4; k++){
			for(int l = 0; l < 4; l++){
				cubicCij(j, k) += cij(i, 16*j + 4*l + k)*xvec[l];
			}
		}
	}

    return CubicInterpolator(y0, dy, ny, cubicCij);
}

CubicInterpolator BicubicInterpolator::reduce_y(const double y){
    int j = findYInterval(y);
    double ybar = (y - y0 - j*dy)/dy;
	Matrix cubicCij(nx, 4);
	double yvec[4] = {1., ybar, ybar*ybar, ybar*ybar*ybar};

	// zj = xi*cij
	for(int i = 0; i < nx; i++){
		for(int k = 0; k < 4; k++){
			for(int l = 0; l < 4; l++){
				cubicCij(i, k) += cij(i, 16*j + 4*k + l)*yvec[l];
			}
		}
	}

    return CubicInterpolator(x0, dx, nx, cubicCij);
}

//////////////////////////////////////////////////////////////////
//////////////           CubicInterpolator        ////////////////
//////////////////////////////////////////////////////////////////

CubicInterpolator::CubicInterpolator(double x0, double dx, const Vector &y): dx(dx), nintervals(y.size()-1), x0(x0), cij(y.size()-1, 4) {
	computeSplineCoefficients(dx, y);
}

CubicInterpolator::CubicInterpolator(const Vector &x, const Vector &y): dx(x[1] - x[0]), nintervals(x.size() - 1), x0(x[0]), cij(x.size() - 1, 4) {
	if(x.size() != y.size()){
		std::cout << "ERROR: Size of x and y vectors do not match \n";
	}
	computeSplineCoefficients(dx, y);
}

CubicInterpolator::CubicInterpolator(double x0, double dx, int nx, Matrix cij): dx(dx), nintervals(nx), x0(x0), cij(cij) {}

double CubicInterpolator::evaluate(const double x){
	int i = findInterval(x);
	return evaluateInterval(i, x);
}

double CubicInterpolator::derivative(const double x){
	int i = findInterval(x);
	return evaluateDerivativeInterval(i, x)/dx;
}

double CubicInterpolator::derivative2(const double x){
	int i = findInterval(x);
	return evaluateSecondDerivativeInterval(i, x)/(dx*dx);
}

void CubicInterpolator::computeSplineCoefficients(double dx, const Vector &y){
	Vector diffopCentral = {1./12., -8./12., 0., 8./12., -1./12.};
	Vector diffopForward = {-25./12., 4., -3., 4./3., -1./4};
	Vector ydx(y.size());

	for(size_t i = 0; i < ydx.size(); i++){
		for(size_t k = 0; k < diffopForward.size(); k++){
			ydx[i] += y[i + k]*diffopForward[k];
		}
	}

	for(size_t i = 2; i < ydx.size(); i++){
		for(size_t k = 0; k < diffopCentral.size(); k++){
			ydx[i] += y[i - 2 + k]*diffopCentral[k];
		}
	}
	
	for(int i = 0; i < nintervals; i++){
		cij(i, 0) = y[i];
		cij(i, 1) = ydx[i];
		cij(i, 2) = -3.*y[i] + 3.*y[i + 1] - 2.*ydx[i] - ydx[i + 1];
		cij(i, 3) = 2.*y[i] - 2.*y[i + 1] + ydx[i] + ydx[i + 1];
	}
}

double CubicInterpolator::evaluateInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	return cij(i, 0) + cij(i, 1)*xbar + cij(i, 2)*xbar*xbar + cij(i, 3)*xbar*xbar*xbar;
}

double CubicInterpolator::evaluateDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	return cij(i, 1) + 2.*cij(i, 2)*xbar + 3.*cij(i, 3)*xbar*xbar;
}

double CubicInterpolator::evaluateSecondDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	return 2.*cij(i, 2) + 6.*cij(i, 3)*xbar;
}

int CubicInterpolator::findInterval(const double x){
	int i = static_cast<int>((x-x0)/dx);
    if(i >= nintervals){
        return nintervals - 1;
    }
	if(i < 0){
		return 0;
	}
	return i;
}
