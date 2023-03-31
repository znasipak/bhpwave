#include "spline.hpp"
#include <iostream>

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

//////////////////////////////////////////////////////////////////
//////////////          BicubicInterpolator       ////////////////
//////////////////////////////////////////////////////////////////

BicubicInterpolator::BicubicInterpolator(const EigenArray &x, const EigenArray &y, const EigenArray &z): BicubicInterpolator(x[0], x[1] - x[0], x.size() - 1, y[0], y[1] - y[0], y.size() - 1, z.matrix()) {}
BicubicInterpolator::BicubicInterpolator(const EigenVector &x, const EigenVector &y, const Eigen::MatrixXd &z): BicubicInterpolator(x[0], x[1] - x[0], x.size() - 1, y[0], y[1] - y[0], y.size() - 1, z) {}
BicubicInterpolator::BicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Eigen::MatrixXd &z): dx(dx), dy(dy), nx(nx), ny(ny), x0(x0), y0(y0), cij(4*nx, 4*ny) {
	if(nx + 1 != z.rows() && ny + 1 != z.cols()){
		if(nx + 1 == z.cols() && ny + 1 == z.rows()){
			// switch x and y
			cij.transposeInPlace();
			computeSplineCoefficients(z);
		}else if((nx + 1)*(ny + 1) == z.size()){
			Eigen::MatrixXd m_z = z.reshaped(ny + 1, nx + 1).transpose();
			computeSplineCoefficients(m_z);
		}else{
			std::cout << "ERROR: Indices of vectors and matrices do not match \n";
		}
	}else{
		computeSplineCoefficients(z);
	}
	std::cout << "Computed spline for z[0, 0] = "<<z(0,0)<<" \n";
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

Eigen::MatrixXd BicubicInterpolator::evaluate(const Eigen::VectorXd &x, const Eigen::VectorXd &y){
    Eigen::MatrixXd mat(x.size(), y.size());

    // this returns all of the intervals that the x and y values are on
	Eigen::ArrayXi iarray = findXInterval(x);
	Eigen::ArrayXi jarray = findYInterval(y);
    
    // at the moment this is slower than just iterating over the x and y lists with nested for loops.
    // I'm not quite sure how to optimize this.

    // first we see how many of the x and y values lie in the same interval
    int jdxi = 0; // index of the first value that lies in the interval
    int jdxf = jdxi; // index of the last value that lies in the interval
    int jnext = jarray(jdxf); // interval of the next y value
    int j = jnext; // interval of the current y value
    while(jdxf < jarray.size() - 1){
        int idxi = 0; // index of the first value that lies in the interval
        int idxf = idxi; // index of the last value that lies in the interval
        int inext = iarray(idxf); // interval of the next x value
        int i = inext; // interval of the current x value
        // std::cout << "Searching with final index ("<<idxf<<", "<<jdxf<<") \n";
        while(j == jnext && jdxf < jarray.size() - 1){// keep searching until the current y value and next y value are on different intervals
            jdxf++;
            jnext = jarray(jdxf);
            while(idxf < iarray.size() - 1){
                while(i == inext && idxf < iarray.size() - 1){ // keep searching until the current x value and next x value are on different intervals
                    idxf++;
                    inext = iarray(idxf);
                }
                // std::cout << "Block starting at ("<<idxi<<", "<<jdxi<<"), x idx in ["<<idxi<<", "<<idxf-1<<"], y idx in ["<<jdxi<<", "<<jdxf-1<<"]\n";
                // std::cout << evaluateInterval(i, j, x(Eigen::seq(idxi, idxf-1)), y(Eigen::seq(jdxi, jdxf-1))) << "\n";
                // std::cout << mat.block(idxi, jdxi, idxf-idxi, jdxf-jdxi) << "\n";
                mat.block(idxi, jdxi, idxf-idxi, jdxf-jdxi) = evaluateInterval(i, j, x(Eigen::seq(idxi, idxf-1)), y(Eigen::seq(jdxi, jdxf-1)));
                // std::cout << mat.block(idxi, jdxi, idxf-idxi, jdxf-jdxi) << "\n";
                idxi = idxf;
                i = inext;
            }
        }
        jdxi = jdxf;
        j = jnext;
    }
	return mat;
}

Eigen::MatrixXd BicubicInterpolator::computeXDerivatives(Eigen::MatrixXd m_z){
    int nsize = m_z.rows();
    Eigen::MatrixXd diffopx = Eigen::MatrixXd::Zero(nsize, nsize);
	diffopx.diagonal(1) = Eigen::VectorXd::Constant(nsize - 1, 8./12.);
	diffopx.diagonal(-1) = Eigen::VectorXd::Constant(nsize - 1, -8./12.);
    diffopx.diagonal(2) = Eigen::VectorXd::Constant(nsize - 2, -1./12.);
	diffopx.diagonal(-2) = Eigen::VectorXd::Constant(nsize - 2, 1./12.);

    // 4th-order forward difference coefficients
    Eigen::Matrix<double, 1, 5> coeffs;
    coeffs << -25./12., 4., -3., 4./3., -1./4;
    for(int i = 0; i < 2; i++){
        diffopx.row(i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(i, Eigen::seq(i, i + 4)) = coeffs;

        diffopx.row(nsize - 1 - i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(nsize - 1 - i, Eigen::seq(nsize - 5 - i, nsize - 1 - i)) = -coeffs.reverse();
    }

	// generate x derivatives
	Eigen::MatrixXd m_zdx = diffopx*m_z/dx;
    return m_zdx;
}

Eigen::MatrixXd BicubicInterpolator::computeYDerivatives(Eigen::MatrixXd m_z){
    int nsize = m_z.cols();
    Eigen::MatrixXd diffopy = Eigen::MatrixXd::Zero(nsize, nsize);
    diffopy.diagonal(1) = Eigen::VectorXd::Constant(nsize - 1, 8./12.);
	diffopy.diagonal(-1) = Eigen::VectorXd::Constant(nsize - 1, -8./12.);
    diffopy.diagonal(2) = Eigen::VectorXd::Constant(nsize - 2, -1./12.);
	diffopy.diagonal(-2) = Eigen::VectorXd::Constant(nsize - 2, 1./12.);

    // 4th-order forward difference coefficients
    Eigen::Matrix<double, 1, 5> coeffs;
    coeffs << -25./12., 4., -3., 4./3., -1./4;
    for(int i = 0; i < 2; i++){
        diffopy.row(i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopy(i, Eigen::seq(i, i + 4)) = coeffs;

        diffopy.row(nsize - 1 - i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopy(nsize - 1 - i, Eigen::seq(nsize - 5 - i, nsize - 1 - i)) = -coeffs.reverse();
    }
    // std::cout << diffopy << "\n";
    diffopy.transposeInPlace();

	// generate y derivatives
	Eigen::MatrixXd m_zdy = m_z*diffopy/dy;
    return m_zdy;
}

Eigen::MatrixXd BicubicInterpolator::computeXYDerivatives(Eigen::MatrixXd m_z){
    Eigen::MatrixXd diffopx = Eigen::MatrixXd::Zero(m_z.rows(), m_z.rows());
	Eigen::MatrixXd diffopy = Eigen::MatrixXd::Zero(m_z.cols(), m_z.cols());
	diffopx.diagonal(1) = Eigen::VectorXd::Constant(m_z.rows() - 1, 0.5);
	diffopx.diagonal(-1) = Eigen::VectorXd::Constant(m_z.rows() - 1, -0.5);
	diffopy.diagonal(-1) = Eigen::VectorXd::Constant(m_z.cols() - 1, 0.5);
	diffopy.diagonal(1) = Eigen::VectorXd::Constant(m_z.cols() - 1, -0.5);
	diffopx(0, 0) = -1.;
	diffopx(0, 1) = 1.;
	diffopx(diffopx.rows() - 1, diffopx.cols() - 2) = -1.;
	diffopx(diffopx.rows() - 1, diffopx.cols() - 1) = 1.;
	diffopy(0, 0) = -1.;
	diffopy(1, 0) = 1.;
	diffopy(diffopy.rows() - 2, diffopy.cols() - 1) = -1.;
	diffopy(diffopy.rows() - 1, diffopy.cols() - 1) = 1.;
    // std::cout << diffopy << "\n";

	// generate y derivatives
	Eigen::MatrixXd m_zdy = m_z*diffopy/dy;
	// generate xy derivatives
	Eigen::MatrixXd m_zdxdy = diffopx*m_zdy/dx;
    return m_zdxdy;
}

Eigen::MatrixXd BicubicInterpolator::computeDerivatives(Eigen::MatrixXd m_z){
	// we approximate first derivatives using central finite differences
	//		f_x(x, y) = [f(x + h, y) - f(x - h, y)]/(2h)
	//		f_y(x, y) = [f(x, y + h) - f(x, y - h)]/(2h)
	// and similarly for the second derivatives
	//		f_xy(x, y) = [f(x + h, y + k) - f(x + h, y - k) - f(x - h, y + k) + f(x - h, y - k)]/(2hk)
	
	// define difference operators in x and y directions
    int nsize = m_z.rows();
    Eigen::MatrixXd diffopx = Eigen::MatrixXd::Zero(nsize, nsize);
	diffopx.diagonal(1) = Eigen::VectorXd::Constant(nsize - 1, 8./12.);
	diffopx.diagonal(-1) = Eigen::VectorXd::Constant(nsize - 1, -8./12.);
    diffopx.diagonal(2) = Eigen::VectorXd::Constant(nsize - 2, -1./12.);
	diffopx.diagonal(-2) = Eigen::VectorXd::Constant(nsize - 2, 1./12.);

    // 4th-order forward difference coefficients
    Eigen::Matrix<double, 1, 5> coeffs;
    coeffs << -25./12., 4., -3., 4./3., -1./4;
    for(int i = 0; i < 2; i++){
        diffopx.row(i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(i, Eigen::seq(i, i + 4)) = coeffs;

        diffopx.row(nsize - 1 - i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(nsize - 1 - i, Eigen::seq(nsize - 5 - i, nsize - 1 - i)) = -coeffs.reverse();
    }
    
    nsize = m_z.cols();
    Eigen::MatrixXd diffopy = Eigen::MatrixXd::Zero(nsize, nsize);
    diffopy.diagonal(1) = Eigen::VectorXd::Constant(nsize - 1, 8./12.);
	diffopy.diagonal(-1) = Eigen::VectorXd::Constant(nsize - 1, -8./12.);
    diffopy.diagonal(2) = Eigen::VectorXd::Constant(nsize - 2, -1./12.);
	diffopy.diagonal(-2) = Eigen::VectorXd::Constant(nsize - 2, 1./12.);

    // 4th-order forward difference coefficients
    for(int i = 0; i < 2; i++){
        diffopy.row(i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopy(i, Eigen::seq(i, i + 4)) = coeffs;

        diffopy.row(nsize - 1 - i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopy(nsize - 1 - i, Eigen::seq(nsize - 5 - i, nsize - 1 - i)) = -coeffs.reverse();
    }
    diffopy.transposeInPlace();

	// generate x derivatives
	Eigen::MatrixXd m_zdx = diffopx*m_z;
	// generate y derivatives
	Eigen::MatrixXd m_zdy = m_z*diffopy;
	// generate xy derivatives
	Eigen::MatrixXd m_zdxdy = diffopx*m_zdy;

    // std::cout << "Computed derivatives\n";
	Eigen::MatrixXd derivmat(cij.rows(), cij.cols());

	for(int j = 0; j < ny; j++){
		for(int i = 0; i < nx; i++){
			Eigen::Matrix4d dmat;
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
            // std::cout << i << ", " << j << "\n";
            // std::cout << dmat << "\n";
			derivmat.block<4,4>(4*i, 4*j) = dmat;
		}
	}

	return derivmat;
}

void BicubicInterpolator::computeSplineCoefficients(Eigen::MatrixXd m_z){
	Eigen::Matrix4d lmat;
	Eigen::Matrix4d rmat;
	Eigen::MatrixXd dmat;

	lmat << 1, 0, 0, 0,
			0, 0, 1, 0,
			-3, 3, -2, -1,
			2, -2, 1, 1;
	rmat << 1, 0, -3, 2,
			0, 0, 3, -2,
			0, 1, -2, 1,
			0, 0, -1, 1;
	dmat = computeDerivatives(m_z);

	for(int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			cij.block<4,4>(4*i, 4*j) = lmat*dmat.block<4,4>(4*i, 4*j)*rmat;
		}
	}
}

double BicubicInterpolator::evaluateInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	Eigen::Vector4d yvec = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

double BicubicInterpolator::evaluateDerivativeXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {0., 1., 2.*xbar, 3.*xbar*xbar};
	Eigen::Vector4d yvec = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

double BicubicInterpolator::evaluateDerivativeYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	Eigen::Vector4d yvec = {0., 1., 2.*ybar, 3.*ybar*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

double BicubicInterpolator::evaluateDerivativeXYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {0., 1., 2.*xbar, 3.*xbar*xbar};
	Eigen::Vector4d yvec = {0., 1., 2.*ybar, 3.*ybar*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

double BicubicInterpolator::evaluateDerivativeXXInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {0., 0., 2., 6.*xbar};
	Eigen::Vector4d yvec = {1, ybar, ybar*ybar, ybar*ybar*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

double BicubicInterpolator::evaluateDerivativeYYInterval(int i, int j, const double x, const double y){
	double xbar = (x - x0 - i*dx)/dx;
	double ybar = (y - y0 - j*dy)/dy;
	Eigen::RowVector4d xvec = {1, xbar, xbar*xbar, xbar*xbar*xbar};
	Eigen::Vector4d yvec = {0., 0., 2., 6.*ybar};
	
	return xvec*(cij.block<4,4>(4*i, 4*j))*yvec;
}

Eigen::MatrixXd BicubicInterpolator::evaluateInterval(int i, int j, const Eigen::VectorXd &x, const Eigen::VectorXd &y){
	Eigen::Array<double, Eigen::Dynamic, 1> xbar = ((x.array() - x0 - i*dx)/dx);
	Eigen::Array<double, 1, Eigen::Dynamic> ybar = ((y.transpose().array() - y0 - j*dy)/dy);
    Eigen::MatrixXd xmat = Eigen::MatrixXd::Constant(xbar.size(), 4, 1.);
    Eigen::MatrixXd ymat = Eigen::MatrixXd::Constant(4, ybar.size(), 1.);
    xmat.col(1) = xbar.matrix();
    xmat.col(2) = xbar.square().matrix();
    xmat.col(3) = xbar.cube().matrix();

    ymat.row(1) = ybar.matrix();
    ymat.row(2) = ybar.square().matrix();
    ymat.row(3) = ybar.cube().matrix();
	
	return xmat*(cij.block<4,4>(4*i, 4*j))*ymat;
}

int BicubicInterpolator::findXInterval(const double x){
	if(x < x0 || x > x0 + nx*dx){
		std::cout << "(ERROR): Value out of bounds \n";
        return 0;
	}
	int i = static_cast<int>((x-x0)/dx);
    if(i == nx){
        return i - 1;
    }
	return i;
}

Eigen::ArrayXi BicubicInterpolator::findXInterval(const Eigen::VectorXd &x){
	if(x(0) < x0 || x(x.size() - 1) >  x0 + nx*dx){
		std::cout << "(ERROR): Value out of bounds \n";
        return Eigen::ArrayXi(x.size());
	}
	Eigen::ArrayXi i = ((x.array() - x0)/dx).floor().cast<int>();
    if(x(x.size() - 1) == x0 + nx*dx){
        i(i.size() - 1) -= 1; // account for the case where we need to evaluate the interpolant at exactly the last boundary point
    }
	return i;
}

int BicubicInterpolator::findYInterval(const double y){
	if(y < y0 || y > y0 + ny*dy){
		std::cout << "(ERROR): Value out of bounds \n";
        return 0;
	}
	int i = static_cast<int>((y-y0)/dy);
    if(i == ny){
        return i - 1;
    }
	return i;
}

Eigen::ArrayXi BicubicInterpolator::findYInterval(const Eigen::VectorXd &y){
	if(y(0) < y0 || y(y.size() - 1) > y0 + ny*dy){
		std::cout << "(ERROR): Value out of bounds \n";
        return Eigen::ArrayXi(y.size());
	}
	Eigen::ArrayXi i = ((y.array()-y0)/dy).floor().cast<int>();
    if(y(y.size() - 1) == y0 + ny*dy){
        i(i.size() - 1) -= 1; // account for the case where we need to evaluate the interpolant at exactly the last boundary point
    }
	return i;
}

CubicInterpolator BicubicInterpolator::reduce_x(const double x){
    int i = findXInterval(x);
    double xbar = (x - x0 - i*dx)/dx;
    Eigen::RowVector4d xvec = {1, xbar, xbar*xbar, xbar*xbar*xbar};

    Eigen::MatrixXd cubicCij(4, ny);
    for(int j = 0; j < ny; j++){
        cubicCij.col(j) = (xvec*cij.block<4,4>(4*i, 4*j)).transpose();
    }

    return CubicInterpolator(dy, cubicCij);
}

CubicInterpolator BicubicInterpolator::reduce_y(const double y){
    int i = findYInterval(y);
    double ybar = (y - y0 - i*dy)/dy;
    Eigen::Vector4d yvec = {1, ybar, ybar*ybar, ybar*ybar*ybar};

    Eigen::MatrixXd cubicCij(4, nx);
    for(int j = 0; j < nx; j++){
        cubicCij.col(j) = (cij.block<4,4>(4*j, 4*i)*yvec);
    }

    return CubicInterpolator(dx, cubicCij);
}

//////////////////////////////////////////////////////////////////
//////////////           CubicInterpolator        ////////////////
//////////////////////////////////////////////////////////////////

CubicInterpolator::CubicInterpolator(double x0, double dx, const EigenArray &y): dx(dx), nintervals(y.size()-1), x0(x0), cij(4, y.size()-1) {
	computeSplineCoefficients(dx, y);
}

CubicInterpolator::CubicInterpolator(const EigenArray &x, const EigenArray &y): dx(x[1] - x[0]), nintervals(x.size() - 1), x0(x[0]), cij(4, x.size() - 1) {
	if(x.size() != y.size()){
		std::cout << "ERROR: Size of x and y vectors do not match \n";
	}
	computeSplineCoefficients(dx, y);
}

CubicInterpolator::CubicInterpolator(double x0, double dx, const Eigen::VectorXd &y): dx(dx), nintervals(y.size()-1), x0(x0), cij(4, y.size()-1) {
	computeSplineCoefficients(dx, y.array());
}

CubicInterpolator::CubicInterpolator(const Eigen::VectorXd &x, const Eigen::VectorXd &y): dx(x[1] - x[0]), nintervals(x.size() - 1), x0(x[0]), cij(4, x.size() - 1) {
	if(x.size() != y.size()){
		std::cout << "ERROR: Size of x and y vectors do not match \n";
	}
	computeSplineCoefficients(dx, y.array());
}

CubicInterpolator::CubicInterpolator(double dx, const Eigen::MatrixXd &cij): dx(dx), cij(cij) {}

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

Eigen::VectorXd CubicInterpolator::computeDerivatives(double dx, const EigenVector &y){
    int nsize = y.size();
    Eigen::MatrixXd diffopx = Eigen::MatrixXd::Zero(nsize, nsize);
	diffopx.diagonal(1) = Eigen::VectorXd::Constant(nsize - 1, 8./12.);
	diffopx.diagonal(-1) = Eigen::VectorXd::Constant(nsize - 1, -8./12.);
    diffopx.diagonal(2) = Eigen::VectorXd::Constant(nsize - 2, -1./12.);
	diffopx.diagonal(-2) = Eigen::VectorXd::Constant(nsize - 2, 1./12.);

    // 4th-order forward difference coefficients
    Eigen::Matrix<double, 1, 5> coeffs;
    coeffs << -25./12., 4., -3., 4./3., -1./4;
    for(int i = 0; i < 2; i++){
        diffopx.row(i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(i, Eigen::seq(i, i + 4)) = coeffs;

        diffopx.row(nsize - 1 - i) = Eigen::MatrixXd::Zero(1, nsize); // zero out row
        diffopx(nsize - 1 - i, Eigen::seq(nsize - 5 - i, nsize - 1 - i)) = -coeffs.reverse();
    }

	// generate x derivatives
	Eigen::VectorXd m_ydx = diffopx*y/dx;
    return m_ydx;
}

Eigen::MatrixXd CubicInterpolator::computeDerivativeVector(double dx, const Eigen::VectorXd &y){
	// generate x derivatives
	Eigen::VectorXd m_ydx = computeDerivatives(dx, y).transpose();

    // std::cout << "Computed derivatives\n";
	Eigen::MatrixXd derivvec(4, m_ydx.size());

	for(int j = 0; j < m_ydx.size() - 1; j++){
        derivvec(0, j) = y(j);
        derivvec(1, j) = y(j + 1);
        derivvec(2, j) = m_ydx(j);
        derivvec(3, j) = m_ydx(j + 1);
	}

	return derivvec;
}

void CubicInterpolator::computeSplineCoefficients(double dx, const EigenArray &y){
	Eigen::Matrix4d lmat;

	lmat << 1, 0, 0, 0,
			0, 0, 0, 1,
			-3, 3, -1, -2,
			2, -2, 1, 1;
	Eigen::MatrixXd dmat = computeDerivativeVector(dx, y.matrix());
    cij = lmat*dmat;
}

double CubicInterpolator::evaluateInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	Eigen::RowVector4d xvec = {1., xbar, xbar*xbar, xbar*xbar*xbar};
	
	return xvec*cij.col(i);
}

double CubicInterpolator::evaluateDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	Eigen::RowVector4d xvec = {0., 1., 2.*xbar, 3.*xbar*xbar};
	
	return xvec*cij.col(i);
}

double CubicInterpolator::evaluateSecondDerivativeInterval(int i, const double x){
	double xbar = (x - x0 - i*dx)/dx;
	Eigen::RowVector4d xvec = {0., 0., 2., 6.*xbar};
	
	return xvec*cij.col(i);
}

int CubicInterpolator::findInterval(const double x){
	if(x < x0 || x > x0 + nintervals*dx){
		std::cout << "(ERROR): Value out of bounds \n";
        return 0;
	}
	int i = static_cast<int>((x-x0)/dx);
    if(i == nintervals){
        return i - 1;
    }
	return i;
}