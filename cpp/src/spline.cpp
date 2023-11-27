#include "spline.hpp"
#include <iostream>

#define ENDPOINT_TOL 1.e-10

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
		_A[i*_m + j] = row[j];
	}
}
void Matrix::col_replace(int i, Vector col){
	for(int j = 0; j < _n; j++){
		_A[j*_m + i] = col[j];
	}
}

Vector Matrix::row(int i){
	Vector row(_m);
	for(int j = 0; j < _m; j++){
		row[j] = _A[i*_m + j];
	}
	return row;
}
Vector Matrix::col(int i){
	Vector col(_n);
	for(int j = 0; j < _n; j++){
		col[j] = _A[j*_m + i];
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

void Matrix::set_value(int i, int j, double val){
	_A[i*_m + j] = val;
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
//////////////          BicubicSpline       ////////////////
//////////////////////////////////////////////////////////////////

BicubicSpline::BicubicSpline(const Vector &x, const Vector &y, Matrix &z, int method): BicubicSpline(x[0], x[1] - x[0], x.size() - 1, y[0], y[1] - y[0], y.size() - 1, z) {}
BicubicSpline::BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, Matrix &z, int method): dx(dx), dy(dy), nx(nx), ny(ny), x0(x0), y0(y0), cij(nx, 16*ny) {
	if(nx + 1 != z.rows() && ny + 1 != z.cols()){
		if(nx + 1 == z.cols() && ny + 1 == z.rows()){
			// switch x and y
			cij.transposeInPlace();
			computeSplineCoefficients(z, method);
		}else if((nx + 1)*(ny + 1) == z.size()){
			Matrix m_z = z.reshaped(ny + 1, nx + 1).transpose();
			computeSplineCoefficients(m_z, method);
		}else{
			std::cout << "ERROR: Indices of vectors and matrices do not match \n";
		}
	}else{
		computeSplineCoefficients(z, method);
	}
}

BicubicSpline::BicubicSpline(const Vector &x, const Vector &y, const Vector &z, int method): BicubicSpline(x[0], x[1] - x[0], x.size() - 1, y[0], y[1] - y[0], y.size() - 1, z) {}
BicubicSpline::BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, const Vector &z_vec, int method): dx(dx), dy(dy), nx(nx), ny(ny), x0(x0), y0(y0), cij(nx, 16*ny) {
	Matrix z(nx+1, ny+1, z_vec);
	if(nx + 1 != z.rows() && ny + 1 != z.cols()){
		if(nx + 1 == z.cols() && ny + 1 == z.rows()){
			// switch x and y
			cij.transposeInPlace();
			computeSplineCoefficients(z, method);
		}else if((nx + 1)*(ny + 1) == z.size()){
			Matrix m_z = z.reshaped(ny + 1, nx + 1).transpose();
			computeSplineCoefficients(m_z, method);
		}else{
			std::cout << "ERROR: Indices of vectors and matrices do not match \n";
		}
	}else{
		computeSplineCoefficients(z, method);
	}
}

double BicubicSpline::evaluate(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateInterval(i, j, x, y);
}

double BicubicSpline::derivative_x(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXInterval(i, j, x, y)/dx;
}

double BicubicSpline::derivative_y(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeYInterval(i, j, x, y)/dy;
}

double BicubicSpline::derivative_xy(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXYInterval(i, j, x, y)/dx/dy;
}

double BicubicSpline::derivative_xx(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeXXInterval(i, j, x, y)/dx/dx;
}

double BicubicSpline::derivative_yy(const double x, const double y){
	int i = findXInterval(x);
	int j = findYInterval(y);
	return evaluateDerivativeYYInterval(i, j, x, y)/dy/dy;
}

Matrix BicubicSpline::computeSplineCoefficientsDY(Matrix &m_z, int method){
	int Nx = m_z.rows();
	int Ny = m_z.cols();
	Matrix m_zdy(Nx, Ny);
	for(int i = 0; i < Nx; i++){
		Vector z_xi = m_z.row(i);
		CubicSpline f_xi = CubicSpline(y0, dy, z_xi, method);
		for(int j = 0; j < Ny; j++){
			m_zdy(i, j) = dy*f_xi.derivative(y0 + j*dy);
		}
	}
	return m_zdy; 
}

Matrix BicubicSpline::computeSplineCoefficientsDX(Matrix &m_z, int method){
	int Nx = m_z.rows();
	int Ny = m_z.cols();
	Matrix m_zdx(Nx, Ny);
	for(int j = 0; j < Ny; j++){
		Vector z_yj = m_z.col(j);
		CubicSpline f_yj = CubicSpline(x0, dx, z_yj, method);
		for(int i = 0; i < Nx; i++){
			m_zdx(i, j) = dx*f_yj.derivative(x0 + i*dx);
		}
	}
	return m_zdx; 
}

void BicubicSpline::computeSplineCoefficients(Matrix &m_z, int method){
	// StopWatch watch;

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
	
	Matrix m_zdx = computeSplineCoefficientsDX(m_z, method);
	Matrix m_zdy = computeSplineCoefficientsDY(m_z, method);
	Matrix m_zdxdy = computeSplineCoefficientsDY(m_zdx, method);
	// Matrix m_zdxdy2 = computeSplineCoefficientsDX(m_zdy);

	// int Nx = m_z.rows();
	// int Ny = m_z.cols();
	// for(int i = 0; i < Nx; i++){
	// 	for(int j = 0; j < Ny; j++){
	// 		// if(j == 0){
	// 		// 	std::cout << "dx = " << m_zdx(i, j) << "\n";
	// 		// }
	// 		// if(i == 0){
	// 		// 	std::cout << "dy = " << m_zdy(i, j) << "\n";
	// 		// }
	// 		std::cout << m_zdxdy2(i, j) << "\n";
	// 		std::cout << m_zdxdy(i, j) << "\n";
	// 	}
	// }

	// now this part we just have to accept as being inefficient because we
	// are mixing rows and columns no matter what. The important thing is that
	// we will store the relevant cofficients close to one another in memory

	// watch.start();
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
								Dmat(k, l) += dmat(k, m)*lmat(l, m); // need transpose of lmat
							}
						}
					}
					for(int k = 0; k < 4; k++){
						for(int l = 0; l < 4; l++){
							for(int m = 0; m < 4; m++){
								cij(i, 16*j + 4*k + l) += lmat(k, m)*Dmat(m, l);
							}
						}
					}
				}
			}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

// void BicubicSpline::computeSplineCoefficients(Matrix &m_z){
// 	StopWatch watch;

// 	watch.start();
// 	Matrix lmat(4, 4, 0.);
// 	lmat(0, 0) = 1.;
// 	lmat(1, 2) = 1.;
// 	lmat(2, 0) = -3.;
// 	lmat(2, 1) = 3.;
// 	lmat(2, 2) = -2.;
// 	lmat(2, 3) = -1.;
// 	lmat(3, 0) = 2.;
// 	lmat(3, 1) = -2.;
// 	lmat(3, 2) = 1.;
// 	lmat(3, 3) = 1.;
	
// 	// 4th-order central difference coefficients
// 	Vector diffopCentral = {1./12., -8./12., 0., 8./12., -1./12.};
// 	// 4th-order forward difference coefficients
// 	Vector diffopForward = {-25./12., 4., -3., 4./3., -1./4};
// 	Matrix m_zdy(m_z.rows(), m_z.cols());

// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdy.rows(); i++){
// 				for(int j = 0; j < 2; j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_z(i, j + k)*diffopForward[k];
// 					}
// 					m_zdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdy.rows(); i++){
// 				for(int j = 2; j < m_zdy.cols() - 2; j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_z(i, j - 2 + k)*diffopCentral[k];
// 					}
// 					m_zdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdy.rows(); i++){
// 				for(int j = m_zdy.cols() - 2; j < m_zdy.cols(); j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += -m_z(i, j - k)*diffopForward[k];
// 					}
// 					m_zdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	watch.stop();
// 	watch.print();
// 	watch.reset();

// 	watch.start();
// 	// Might be more efficient to transpose m_z here to take advantage of
// 	// row-major ordering
// 	Matrix m_zT = m_z.transpose();
// 	Matrix m_zdxT(m_zT.rows(), m_zT.cols());

// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 		for(int i = 0; i < m_zdxT.rows(); i++){
// 			for(int j = 0; j < 2; j++){
// 				double sum = 0.;
// 				for(int k = 0; k < 5; k++){
// 					sum += m_zT(i, j + k)*diffopForward[k];
// 				}
// 				m_zdxT(i, j) = sum;
// 			}
// 		}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdxT.rows(); i++){
// 				for(int j = 2; j < m_zdxT.cols() - 2; j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_zT(i, j - 2 + k)*diffopCentral[k];
// 					}
// 					m_zdxT(i, j) = sum;
// 				}
// 			}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdxT.rows(); i++){
// 				for(int j = m_zdxT.cols() - 2; j < m_zdxT.cols(); j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += -m_zT(i, j - k)*diffopForward[k];
// 					}
// 					m_zdxT(i, j) = sum;
// 				}
// 			}
// 	}
// 	Matrix m_zdx = m_zdxT.transpose();
// 	watch.stop();
// 	watch.print();
// 	watch.reset();
	
// 	watch.start();
// 	Matrix m_zdxdy(m_zdx.rows(), m_zdx.cols());
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdxdy.rows(); i++){
// 				for(int j = 0; j < 2; j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_zdx(i, j + k)*diffopForward[k];
// 					}
// 					m_zdxdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdxdy.rows(); i++){
// 				for(int j = 2; j < m_zdxdy.cols() - 2; j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_zdx(i, j - 2 + k)*diffopCentral[k];
// 					}
// 					m_zdxdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 			for(int i = 0; i < m_zdxdy.rows(); i++){
// 				for(int j = m_zdxdy.cols() - 2; j < m_zdxdy.cols(); j++){
// 					double sum = 0.;
// 					for(int k = 0; k < 5; k++){
// 						sum += m_zdx(i, j - k)*diffopForward[k];
// 					}
// 					m_zdxdy(i, j) = sum;
// 				}
// 			}
// 	}
// 	watch.stop();
// 	watch.print();
// 	watch.reset();

// 	// now this part we just have to accept as being inefficient because we
// 	// are mixing rows and columns no matter what. The important thing is that
// 	// we will store the relevant cofficients close to one another in memory

// 	watch.start();
// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2) schedule(dynamic, 16)
// 			for(int i = 0; i < nx; i++){
// 				for(int j = 0; j < ny; j++){
// 					Matrix dmat(4, 4);
// 					dmat(0, 0) = m_z(i, j); // f(0,0)
// 					dmat(0, 1) = m_z(i, j + 1); // f(0,1)
// 					dmat(0, 2) = m_zdy(i, j); // fy(0,0)
// 					dmat(0, 3) = m_zdy(i, j + 1); // fy(0,1)
// 					dmat(1, 0) = m_z(i + 1, j); // f(1,0)
// 					dmat(1, 1) = m_z(i + 1, j + 1); // f(1,1)
// 					dmat(1, 2) = m_zdy(i + 1, j); // fy(1,0)
// 					dmat(1, 3) = m_zdy(i + 1, j + 1); // fy(1,1)
// 					dmat(2, 0) = m_zdx(i, j); // fx(0,0)
// 					dmat(2, 1) = m_zdx(i, j + 1); // fx(0,1)
// 					dmat(2, 2) = m_zdxdy(i, j); // fxy(0,0)
// 					dmat(2, 3) = m_zdxdy(i, j + 1); // fxy(0,1)
// 					dmat(3, 0) = m_zdx(i + 1, j); // fx(1,0)
// 					dmat(3, 1) = m_zdx(i + 1, j + 1); // fx(1,1)
// 					dmat(3, 2) = m_zdxdy(i + 1, j); // fxy(1,0)
// 					dmat(3, 3) = m_zdxdy(i + 1, j + 1); // fxy(1,1)

// 					// this part is slow. Just lots of matrix multiplication
// 					Matrix Dmat(4, 4);
// 					for(int k = 0; k < 4; k++){
// 						for(int l = 0; l < 4; l++){
// 							for(int m = 0; m < 4; m++){
// 								Dmat(l, k) += dmat(k, m)*lmat(k, m);
// 							}
// 							for(int m = 0; m < 4; m++){
// 								cij(i, 16*j + 4*k + l) += lmat(k, m)*Dmat(k, m);
// 							}
// 						}
// 					}
// 					// for(int k = 0; k < 4; k++){
// 					// 	for(int l = 0; l < 4; l++){
// 					// 		for(int m = 0; m < 4; m++){
// 					// 			cij(i, 16*j + 4*k + l) += lmat(k, m)*Dmat(k, m);
// 					// 		}
// 					// 	}
// 					// }
// 				}
// 			}
// 	}
// 	watch.stop();
// 	watch.print();
// 	watch.reset();
// }

double BicubicSpline::evaluateInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	// double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	// double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	// double zvec[4] = {0., 0., 0., 0.};
	// double result = 0.;
	double zvec[4];
	double result;

	// // zj = cij*yi
	// for(int k = 0; k < 4; k++){
	// 	for(int l = 0; l < 4; l++){
	// 		zvec[k] += cij(i, 16*j + 4*k + l)*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 4; k++){
	// 	result += xvec[k]*zvec[k];
	// }

	for(int k = 0; k < 4; k++){
		zvec[k] = cij(i, 16*j + 4*k + 0) + ybar*(cij(i, 16*j + 4*k + 1) + ybar*(cij(i, 16*j + 4*k + 2) + cij(i, 16*j + 4*k + 3)*ybar));
	}

	result = zvec[0] + xbar*(zvec[1] + xbar*(zvec[2] + zvec[3]*xbar));

	return result;
}

double BicubicSpline::evaluateDerivativeXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	// double xvec[3] = {1., 2.*xbar, 3.*xbar*xbar};
	// double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	// double zvec[3] = {0., 0., 0.};
	// double result = 0.;

	double zvec[4];
	double result;

	// // zj = cij*yi
	// for(int k = 0; k < 3; k++){
	// 	for(int l = 0; l < 4; l++){
	// 		zvec[k] += cij(i, 16*j + 4*(1 + k) + l)*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 3; k++){
	// 	result += xvec[k]*zvec[k];
	// }

	for(int k = 0; k < 4; k++){
		zvec[k] = cij(i, 16*j + 4*k + 0) + ybar*(cij(i, 16*j + 4*k + 1) + ybar*(cij(i, 16*j + 4*k + 2) + cij(i, 16*j + 4*k + 3)*ybar));
	}

	result = (zvec[1] + xbar*(2.*zvec[2] + 3.*zvec[3]*xbar));

	return result;
}

double BicubicSpline::evaluateDerivativeYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;

	// double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	// double yvec[3] = {1., 2.*ybar, 3.*ybar*ybar};
	// double zvec[4] = {0., 0., 0., 0.};
	// double result = 0.;

	double zvec[4];
	double result;

	// zj = cij*yi
	// for(int k = 0; k < 4; k++){
	// 	for(int l = 0; l < 3; l++){
	// 		zvec[k] += cij(i, 16*j + 4*k + (l + 1))*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 4; k++){
	// 	result += xvec[k]*zvec[k];
	// }
	for(int k = 0; k < 4; k++){
		zvec[k] = cij(i, 16*j + 4*k + 1) + ybar*(2.*cij(i, 16*j + 4*k + 2) + 3.*cij(i, 16*j + 4*k + 3)*ybar);
	}

	result = zvec[0] + xbar*(zvec[1] + xbar*(zvec[2] + zvec[3]*xbar));
	
	return result;
}

double BicubicSpline::evaluateDerivativeXYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	// double xvec[3] = {1., 2.*xbar, 3.*xbar*xbar};
	// double yvec[3] = {1., 2.*ybar, 3.*ybar*ybar};
	// double zvec[3] = {0., 0., 0.};
	// double result = 0.;
	double zvec[4];
	double result;

	// zj = cij*yi
	// for(int k = 0; k < 3; k++){
	// 	for(int l = 0; l < 3; l++){
	// 		zvec[k] += cij(i, 16*j + 4*(k + 1) + (l + 1))*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 3; k++){
	// 	result += xvec[k]*zvec[k];
	// }

	for(int k = 0; k < 4; k++){
		zvec[k] = (cij(i, 16*j + 4*k + 1) + ybar*(2.*cij(i, 16*j + 4*k + 2) + 3.*cij(i, 16*j + 4*k + 3)*ybar));
	}

	result = (zvec[1] + xbar*(2.*zvec[2] + 3.*zvec[3]*xbar));
	
	return result;
}

double BicubicSpline::evaluateDerivativeXXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	// double xvec[2] = {2., 6.*xbar};
	// double yvec[4] = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	// double zvec[2] = {0., 0.};
	// double result = 0.;
	double zvec[4];
	double result;

	// // zj = cij*yi
	// for(int k = 0; k < 2; k++){
	// 	for(int l = 0; l < 4; l++){
	// 		zvec[k] += cij(i, 16*j + 4*(2 + k) + l)*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 2; k++){
	// 	result += xvec[k]*zvec[k];
	// }

	for(int k = 0; k < 4; k++){
		zvec[k] = cij(i, 16*j + 4*k + 0) + ybar*(cij(i, 16*j + 4*k + 1) + ybar*(cij(i, 16*j + 4*k + 2) + cij(i, 16*j + 4*k + 3)*ybar));
	}

	result = 2.*(zvec[2] + 3.*zvec[3]*xbar);

	return result;
}

double BicubicSpline::evaluateDerivativeYYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	
	// double xvec[4] = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	// double yvec[2] = {2., 6.*ybar};
	// double zvec[4] = {0., 0., 0., 0.};
	// double result = 0.;
	double zvec[4];
	double result;

	// zj = cij*yi
	// for(int k = 0; k < 4; k++){
	// 	for(int l = 0; l < 2; l++){
	// 		zvec[k] += cij(i, 16*j + 4*k + (l + 2))*yvec[l];
	// 	}
	// }

	// // result = xj*zj
	// for(int k = 0; k < 4; k++){
	// 	result += xvec[k]*zvec[k];
	// }

	for(int k = 0; k < 4; k++){
		zvec[k] = 2.*(cij(i, 16*j + 4*k + 2) + 3.*cij(i, 16*j + 4*k + 3)*ybar);
	}

	result = zvec[0] + xbar*(zvec[1] + xbar*(zvec[2] + zvec[3]*xbar));
	
	return result;
}

int BicubicSpline::findXInterval(const double x){
	int i = static_cast<int>((x-x0)/dx);
    if(i >= nx){
        return nx - 1;
    }
	if( i < 0){
		return 0;
	}
	return i;
}

int BicubicSpline::findYInterval(const double y){
	int i = static_cast<int>((y-y0)/dy);
    if(i >= ny){
        return ny - 1;
    }
	if( i < 0){
		return 0;
	}
	return i;
}

CubicSpline BicubicSpline::reduce_x(const double x){
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

    return CubicSpline(y0, dy, ny, cubicCij);
}

CubicSpline BicubicSpline::reduce_y(const double y){
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

    return CubicSpline(x0, dx, nx, cubicCij);
}

//////////////////////////////////////////////////////////////////
//////////////           CubicSpline        ////////////////
//////////////////////////////////////////////////////////////////

CubicSpline::CubicSpline(double x0, double dx, const Vector &y, int method): dx(dx), nintervals(y.size()-1), x0(x0), cij(y.size()-1, 4) {
	if(method == 0){
		computeSplineCoefficients(dx, y);
	}else if(method == 1){
		computeSplineCoefficientsNotAKnot(dx, y);
	}else if(method == 2){
		computeSplineCoefficientsZeroClamped(dx, y);
	}else if(method == 3){
		computeSplineCoefficientsE3(dx, y);
	}else if(method == 4){
		computeSplineCoefficientsNaturalFirst(dx, y);
	}else{
		computeSplineCoefficientsNotAKnot(dx, y);
	}
}

CubicSpline::CubicSpline(const Vector &x, const Vector &y, int method): dx(x[1] - x[0]), nintervals(x.size() - 1), x0(x[0]), cij(x.size() - 1, 4) {
	if(x.size() != y.size()){
		std::cout << "ERROR: Size of x and y vectors do not match \n";
	}
	if(method == 0){
		computeSplineCoefficients(dx, y);
	}else if(method == 1){
		computeSplineCoefficientsNotAKnot(dx, y);
	}else if(method == 2){
		computeSplineCoefficientsZeroClamped(dx, y);
	}else if(method == 3){
		computeSplineCoefficientsE3(dx, y);
	}else if(method == 4){
		computeSplineCoefficientsNaturalFirst(dx, y);
	}else{
		computeSplineCoefficientsNotAKnot(dx, y);
	}
}

CubicSpline::CubicSpline(double x0, double dx, int nx, Matrix cij): dx(dx), nintervals(nx), x0(x0), cij(cij) {}

double CubicSpline::getSplineCoefficient(int i, int j){
	return cij(i, j);
}

double CubicSpline::evaluate(const double x){
	int i = findInterval(x);
	return evaluateInterval(i, x);
}

double CubicSpline::derivative(const double x){
	int i = findInterval(x);
	return evaluateDerivativeInterval(i, x);
}

double CubicSpline::derivative2(const double x){
	int i = findInterval(x);
	return evaluateSecondDerivativeInterval(i, x);
}

void CubicSpline::computeSplineCoefficients(double dx, const Vector &y){
	// Calculation with natural boundary conditions that follows the GSL algorithm
	// Essentially we first calculate the second-derivatives assuming y''(x0) = y''(xn) = 0
	// Then from the values of y and y'', we construct the spline coefficients

	int nsize = y.size() - 2;
	// Forward sweep method for tridiagonal solve on Wikipedia
	double a = 1.;
	Vector b(nsize, 4.);
	double c = 1.;
	Vector d(nsize);
	Vector ydx2Over2(nsize);

	for(int i = 1; i <= nsize; i++){
		d[i - 1] = 3.*(y[i + 1] - 2.0*y[i] + y[i - 1])/(dx*dx);
	}

	//forward sweep
	double w = 0.;
	for(int i = 2; i <= nsize; i++){
		w = a/b[i - 2];
		b[i - 1] = b[i - 1] - w*c;
		d[i - 1] = d[i - 1] - w*d[i - 2];
	}
	// back substitution
	ydx2Over2[nsize - 1] = d[nsize - 1]/b[nsize - 1];
	for(int i = 1; i < nsize; i++){
		ydx2Over2[nsize - i - 1] = (d[nsize - i - 1] - c*ydx2Over2[nsize - i])/b[nsize - i - 1];
	}

	int i = 0;
	cij(i, 0) = y[i];
	cij(i, 2) = 0;
	cij(i, 1) = (y[i + 1] - y[i])/dx - dx*(ydx2Over2[i])/3.0;
	cij(i, 3) = (ydx2Over2[i])/(3.0*dx);
	
	for(i = 1; i < nintervals - 1; i++){
		cij(i, 0) = y[i];
		cij(i, 2) = ydx2Over2[i - 1];
		cij(i, 1) = (y[i + 1] - y[i])/dx - dx*(ydx2Over2[i] + 2.*ydx2Over2[i - 1])/3.0;
		cij(i, 3) = (ydx2Over2[i] - ydx2Over2[i - 1])/(3.0*dx);
	}

	i = nintervals - 1;
	cij(i, 0) = y[i];
	cij(i, 2) = ydx2Over2[i - 1];
	cij(i, 1) = (y[i + 1] - y[i])/dx - dx*(2.*ydx2Over2[i - 1])/3.0;
	cij(i, 3) = -(ydx2Over2[i - 1])/(3.0*dx);
}

void CubicSpline::computeSplineCoefficientsNaturalFirst(double dx, const Vector &y){
	// Calculation with not-a-know boundary conditions that follows the GSL algorithm
	// Essentially we calculate the first-derivatives assuming y'''(x0) = y'''(x1) and y'''(xn) = y'''(xn-1)
	// Then from the values of y and y', we construct the spline coefficients

	int nsize = y.size();
	// Forward sweep method for tridiagonal solve on Wikipedia
	Vector a(nsize, 1.);
	Vector b(nsize, 4.);
	Vector c(nsize, 1.);
	Vector d(nsize);
	Vector ydx(nsize);

	b[0] = 2.;
	c[0] = 1.;  
	b[nsize - 1] = 2.;
	a[nsize - 1] = 1.;

	for(int i = 1; i < nsize - 1; i++){
		d[i] = 3.*(y[i + 1] - y[i - 1])/(dx);
	}
	d[0] = 3.*(y[1] - y[0])/dx;
	d[nsize - 1] = 3.*(y[nsize - 1] - y[nsize - 2])/dx;

	//forward sweep
	double w = 0.;
	double temp = 0.;
	for(int i = 1; i < nsize; i++){
		w = a[i]/b[i - 1];
		temp = b[i] - w*c[i - 1];
		b[i] = temp;
		temp = d[i] - w*d[i - 1];
		d[i] = temp;
	}
	// back substitution
	ydx[nsize - 1] = d[nsize - 1]/b[nsize - 1];
	for(int i = 1; i < nsize; i++){
		ydx[nsize - i - 1] = (d[nsize - i - 1] - c[nsize - i - 1]*ydx[nsize - i])/b[nsize - i - 1];
	}
	
	for(int i = 0; i < nintervals; i++){
		cij(i, 0) = y[i];
		cij(i, 1) = ydx[i];
		cij(i, 2) = 3.*(y[i + 1] - y[i])/(dx*dx) - (ydx[i + 1] + 2.*ydx[i])/dx;
		cij(i, 3) = (ydx[i + 1] - ydx[i])/(3.*dx*dx) - 2.*cij(i, 2)/(3.*dx);
	}
}

void CubicSpline::computeSplineCoefficientsZeroClamped(double dx, const Vector &y){
	// Calculation with not-a-know boundary conditions that follows the GSL algorithm
	// Essentially we calculate the first-derivatives assuming y'''(x0) = y'''(x1) and y'''(xn) = y'''(xn-1)
	// Then from the values of y and y', we construct the spline coefficients

	int nsize = y.size();
	// Forward sweep method for tridiagonal solve on Wikipedia
	Vector a(nsize, 1.);
	Vector b(nsize, 4.);
	Vector c(nsize, 1.);
	Vector d(nsize);
	Vector ydx(nsize);

	b[0] = 1.;
	c[0] = 0.;  
	b[nsize - 1] = 1.;
	a[nsize - 1] = 0.;

	for(int i = 1; i < nsize - 1; i++){
		d[i] = 3.*(y[i + 1] - y[i - 1])/(dx);
	}
	d[0] = 0.;
	d[nsize - 1] = 0.;

	//forward sweep
	double w = 0.;
	for(int i = 1; i < nsize; i++){
		w = a[i]/b[i - 1];
		b[i] = b[i] - w*c[i - 1];
		d[i] = d[i] - w*d[i - 1];
	}
	// back substitution
	ydx[nsize - 1] = d[nsize - 1]/b[nsize - 1];
	for(int i = 1; i < nsize; i++){
		ydx[nsize - i - 1] = (d[nsize - i - 1] - c[nsize - i - 1]*ydx[nsize - i])/b[nsize - i - 1];
	}
	
	for(int i = 0; i < nintervals; i++){
		cij(i, 0) = y[i];
		cij(i, 1) = ydx[i];
		cij(i, 2) = 3.*(y[i + 1] - y[i])/(dx*dx) - (ydx[i + 1] + 2.*ydx[i])/dx;
		cij(i, 3) = (ydx[i + 1] - ydx[i])/(3.*dx*dx) - 2.*cij(i, 2)/(3.*dx);
	}
}

void CubicSpline::computeSplineCoefficientsE3(double dx, const Vector &y){
	// Calculation with not-a-know boundary conditions that follows the GSL algorithm
	// Essentially we calculate the first-derivatives assuming y'''(x0) = y'''(x1) and y'''(xn) = y'''(xn-1)
	// Then from the values of y and y', we construct the spline coefficients

	int nsize = y.size();
	// Forward sweep method for tridiagonal solve on Wikipedia
	Vector a(nsize, 1.);
	Vector b(nsize, 4.);
	Vector c(nsize, 1.);
	Vector d(nsize);
	Vector ydx(nsize);

	b[0] = 6.;
	c[0] = 18.;  
	b[nsize - 1] = 6.;
	a[nsize - 1] = 18.;

	for(int i = 1; i < nsize - 1; i++){
		d[i] = 3.*(y[i + 1] - y[i - 1])/(dx);
	}
	d[0] = (-y[3] + 9.*y[2] + 9.*y[1] - 17.*y[0])/dx;
	d[nsize - 1] = (y[nsize - 1 - 3] - 9.*y[nsize - 1 - 2] - 9.*y[nsize - 1 - 1] + 17.*y[nsize - 1])/dx;

	//forward sweep
	double w = 0.;
	for(int i = 1; i < nsize; i++){
		w = a[i]/b[i - 1];
		b[i] = b[i] - w*c[i - 1];
		d[i] = d[i] - w*d[i - 1];
	}
	// back substitution
	ydx[nsize - 1] = d[nsize - 1]/b[nsize - 1];
	for(int i = 1; i < nsize; i++){
		ydx[nsize - i - 1] = (d[nsize - i - 1] - c[nsize - i - 1]*ydx[nsize - i])/b[nsize - i - 1];
	}
	
	for(int i = 0; i < nintervals; i++){
		cij(i, 0) = y[i];
		cij(i, 1) = ydx[i];
		cij(i, 2) = 3.*(y[i + 1] - y[i])/(dx*dx) - (ydx[i + 1] + 2.*ydx[i])/dx;
		cij(i, 3) = (ydx[i + 1] - ydx[i])/(3.*dx*dx) - 2.*cij(i, 2)/(3.*dx);
	}
}

void CubicSpline::computeSplineCoefficientsNotAKnot(double dx, const Vector &y){
	// Calculation with not-a-know boundary conditions that follows the GSL algorithm
	// Essentially we calculate the first-derivatives assuming y'''(x0) = y'''(x1) and y'''(xn) = y'''(xn-1)
	// Then from the values of y and y', we construct the spline coefficients

	int nsize = y.size();
	// Forward sweep method for tridiagonal solve on Wikipedia
	Vector a(nsize, 1.);
	Vector b(nsize, 4.);
	Vector c(nsize, 1.);
	Vector d(nsize);
	Vector ydx(nsize);

	b[0] = 2.;
	c[0] = 4.;  
	b[nsize - 1] = 2.;
	a[nsize - 1] = 4.;

	for(int i = 1; i < nsize - 1; i++){
		d[i] = 3.*(y[i + 1] - y[i - 1])/(dx);
	}
	// d[0] = 3.*(y[1] - y[0])/dx;
	// d[nsize - 1] = 3.*(y[nsize - 1] - y[nsize - 2])/dx;
	d[0] = (y[2] + 4.*y[1] - 5.*y[0])/dx;
	d[nsize - 1] = (5.*y[nsize - 1] - 4.*y[nsize - 2] - y[nsize - 3])/dx;

	//forward sweep
	double w = 0.;
	for(int i = 1; i < nsize; i++){
		w = a[i]/b[i - 1];
		b[i] = b[i] - w*c[i - 1];
		d[i] = d[i] - w*d[i - 1];
	}
	// back substitution
	ydx[nsize - 1] = d[nsize - 1]/b[nsize - 1];
	for(int i = 1; i < nsize; i++){
		ydx[nsize - i - 1] = (d[nsize - i - 1] - c[nsize - i - 1]*ydx[nsize - i])/b[nsize - i - 1];
	}
	
	for(int i = 0; i < nintervals; i++){
		cij(i, 0) = y[i];
		cij(i, 1) = ydx[i];
		cij(i, 2) = 3.*(y[i + 1] - y[i])/(dx*dx) - (ydx[i + 1] + 2.*ydx[i])/dx;
		cij(i, 3) = (ydx[i + 1] - ydx[i])/(3.*dx*dx) - 2.*cij(i, 2)/(3.*dx);
	}
}

double CubicSpline::evaluateInterval(int i, const double x){
	double xbar = (x - x0 - i*dx);
	return cij(i, 0) + xbar*(cij(i, 1) + xbar*(cij(i, 2) + cij(i, 3)*xbar));
}

double CubicSpline::evaluateDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx);
	return cij(i, 1) + xbar*(2.0*cij(i, 2) + 3.0*cij(i, 3)*xbar);
}

double CubicSpline::evaluateSecondDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx);
	return 2.0*(cij(i, 2) + 3.0*cij(i, 3)*xbar);
}

int CubicSpline::findInterval(const double x){
	int i = static_cast<int>((x-x0)/dx);
    if(i >= nintervals){
        return nintervals - 1;
    }
	if(i < 0){
		return 0;
	}
	return i;
}
