import numpy as np
cimport numpy as np
from cython.operator import dereference
import os
import warnings
from libcpp.vector cimport vector

cdef extern from "spline.hpp":
    cdef cppclass Matrix:
        Matrix()
        Matrix(int n)
        Matrix(int n, int m)
        Matrix(int n, int m, vector[double] A)
        Matrix(int n, int m, double val)

        void set_value(int i, int j, double val)
        double& operator()(int i, int j)

    cdef cppclass ThreeTensor:
        ThreeTensor()
        ThreeTensor(int nx)
        ThreeTensor(int nx, int ny, int nz)
        ThreeTensor(int nx, int ny, int nz, double *A)
        ThreeTensor(int nx, int ny, int nz, vector[double] A)
        ThreeTensor(int nx, int ny, int nz, double val)

        int rows() const
        int cols() const
        int slcs() const
        int size() const

        void row_replace(int i, Matrix row)
        void col_replace(int j, Matrix col)
        void slc_replace(int k, Matrix slc)

        Matrix row(int i)
        vector[double] rowcol(int i, int j)
        vector[double] rowslc(int i, int k)
        Matrix col(int j)
        vector[double] colslc(int j, int k)
        Matrix slc(int k)

        void reshape(int nx, int ny, int nz)
        ThreeTensor reshaped(int nx, int ny, int nz) const

        void set_value(int i, int j, int k, double val)

        double& operator()(int i, int j, int k)

        vector[double] data()

    cdef cppclass CubicSpline:
        CubicSpline(double x0, double dx, const vector[double] &y, int method) except +
        CubicSpline(const vector[double] &x, const vector[double] &y, int method) except +

        double getSplineCoefficient(int i, int j)

        double evaluate(const double x)
        double derivative(const double x)
        double derivative2(const double x)

    cdef cppclass BicubicSpline:
        BicubicSpline(const vector[double] &x, const vector[double] &y, const Matrix &z, int method) except +
        BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, const Matrix &z, int method) except +
        
        double evaluate(const double x, const double y)
        double derivative_x(const double x, const double y)
        double derivative_y(const double x, const double y)
        double derivative_xy(const double x, const double y)
        double derivative_xx(const double x, const double y)
        double derivative_yy(const double x, const double y)
        CubicSpline reduce_x(const double x)
        CubicSpline reduce_y(const double y)
        double getSplineCoefficient(int i, int j, int nx, int ny)

    cdef cppclass TricubicSpline:
        # TricubicSpline(const vector[double] &x, const vector[double] &y, const vector[double] &z, ThreeTensor &f, int method) except +
        TricubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, double z0, double dz, int nz, ThreeTensor &f, int method) except +
        double evaluate(const double x, const double y, const double z)
        double derivative_x(const double x, const double y, const double z)
        double derivative_y(const double x, const double y, const double z)
        double derivative_z(const double x, const double y, const double z)
        double derivative_xy(const double x, const double y, const double z)
        double derivative_xz(const double x, const double y, const double z)
        double derivative_yz(const double x, const double y, const double z)
        double derivative_xx(const double x, const double y, const double z)
        double derivative_yy(const double x, const double y, const double z)
        double derivative_zz(const double x, const double y, const double z)

        double getSplineCoefficient(int i, int j, int k, int nx, int ny, int nz)
        ThreeTensor getSplineCoefficients()

cdef class CyCubicSpline:
    cdef CubicSpline *scpp

    def __init__(self, double x0, double dx, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] f, int method):
        cdef vector[double] fvec = vector[double](len(f))
        for i in range(len(f)):
            fvec[i] = f[i]
        self.scpp = new CubicSpline(x0, dx, fvec, method)

    def coefficient(self, int i, int j):
        return self.scpp.getSplineCoefficient(i, j)

    def eval(self, double x):
        return self.scpp.evaluate(x)

    def deriv(self, double x):
        return self.scpp.derivative(x)

    def deriv2(self, double x):
        return self.scpp.derivative2(x)
    

cdef class CyBicubicSpline:
    cdef BicubicSpline *scpp

    def __init__(self, double x0, double dx, int nx, double y0, double dy, int ny, np.ndarray[ndim=2, dtype=np.float64_t, mode='c'] f, int method):
        cdef Matrix mz = Matrix(nx + 1, ny + 1)
        for i in range(nx + 1):
            for j in range(ny + 1):
                mz.set_value(i, j, f[i, j])
        self.scpp = new BicubicSpline(x0, dx, nx, y0, dy, ny, mz, method)

    def coefficient(self, int i, int j, int nx, int ny):
        return self.scpp.getSplineCoefficient(i, j, nx, ny)
    
    def eval(self, double x, double y):
        return self.scpp.evaluate(x, y)

    def deriv_x(self, double x, double y):
        return self.scpp.derivative_x(x, y)

    def deriv_y(self, double x, double y):
        return self.scpp.derivative_y(x, y)
    
    def deriv_xx(self, double x, double y):
        return self.scpp.derivative_xx(x, y)

    def deriv_yy(self, double x, double y):
        return self.scpp.derivative_yy(x, y)

    def deriv_xy(self, double x, double y):
        return self.scpp.derivative_xy(x, y)

cdef class CyTricubicSpline:
    cdef TricubicSpline *scpp
    cdef int nx
    cdef int ny
    cdef int nz

    def __init__(self, double x0, double dx, int nx, double y0, double dy, int ny, double z0, double dz, int nz, np.ndarray[ndim=3, dtype=np.float64_t, mode='c'] f, int method):
        cdef int fnx, fny, fnz
        fnx = f.shape[0]
        fny = f.shape[1]
        fnz = f.shape[2]
        self.nx = nx
        self.ny = ny
        self.nz = nz
        cdef ThreeTensor ftens = ThreeTensor(fnx, fny, fnz, &f[0,0,0])
        self.scpp = new TricubicSpline(x0, dx, nx, y0, dy, ny, z0, dz, nz, ftens, method)

    def coefficient(self, int i, int j, int k, int nx, int ny, int nz):
        return self.scpp.getSplineCoefficient(i, j, k, nx, ny, nz)

    def coefficients(self):
        cdef np.ndarray[ndim=3, dtype=np.float64_t, mode='c'] cijk = np.asarray(self.scpp.getSplineCoefficients().data()).reshape(self.nx, self.ny, 64*self.nz)
        return cijk

    def eval(self, double x, double y, double z):
        return self.scpp.evaluate(x, y, z)

    def deriv_x(self, double x, double y, double z):
        return self.scpp.derivative_x(x, y, z)

    def deriv_y(self, double x, double y, double z):
        return self.scpp.derivative_y(x, y, z)

    def deriv_z(self, double x, double y, double z):
        return self.scpp.derivative_z(x, y, z)
    
    def deriv_xx(self, double x, double y, double z):
        return self.scpp.derivative_xx(x, y, z)

    def deriv_yy(self, double x, double y, double z):
        return self.scpp.derivative_yy(x, y, z)

    def deriv_zz(self, double x, double y, double z):
        return self.scpp.derivative_yy(x, y, z)

    def deriv_xy(self, double x, double y, double z):
        return self.scpp.derivative_xy(x, y, z)

    def deriv_xz(self, double x, double y, double z):
        return self.scpp.derivative_xz(x, y, z)

    def deriv_yz(self, double x, double y, double z):
        return self.scpp.derivative_yz(x, y, z)

