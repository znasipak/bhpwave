import numpy as np
cimport numpy as np
from cython.operator import dereference
import os
import warnings
from libcpp.vector cimport vector

cdef extern from "spline.hpp":
    cdef cppclass Spline:
        Spline(vector[double] x, vector[double] f)
        Spline(const Spline& spline)

        double evaluate(const double &x0)
        double derivative(const double &x0)
        double derivative2(const double &x0)

        void reconstruct(vector[double] x, vector[double] f)

    cdef cppclass Spline2D:
        Spline2D(vector[double] x, vector[double] y, vector[double] f)
        Spline2D(const Spline2D& spline)

        double evaluate(const double &x0, const double &y0)
        double derivative_x(const double &x0, const double &y0)
        double derivative_y(const double &x0, const double &y0)
        double derivative_xx(const double &x0, const double &y0)
        double derivative_yy(const double &x0, const double &y0)
        double derivative_xy(const double &x0, const double &y0)

        void reconstruct(vector[double] x, vector[double] y, vector[double] f)

    cdef cppclass Matrix:
        Matrix()
        Matrix(int n)
        Matrix(int n, int m)
        Matrix(int n, int m, vector[double] A)
        Matrix(int n, int m, double val)

        void set_value(int i, int j, double val)
        double& operator()(int i, int j)

    cdef cppclass CubicInterpolator:
        CubicInterpolator(double x0, double dx, const vector[double] &y) except +
        CubicInterpolator(const vector[double] &x, const vector[double] &y) except +

        double evaluate(const double x)
        double derivative(const double x)
        double derivative2(const double x)

    cdef cppclass BicubicInterpolator:
        BicubicInterpolator(const vector[double] &x, const vector[double] &y, const Matrix &z)
        BicubicInterpolator(double x0, double dx, int nx, double y0, double dy, int ny, const Matrix &z)
        
        double evaluate(const double x, const double y)
        double derivative_x(const double x, const double y)
        double derivative_y(const double x, const double y)
        double derivative_xy(const double x, const double y)
        double derivative_xx(const double x, const double y)
        double derivative_yy(const double x, const double y)
        CubicInterpolator reduce_x(const double x)
        CubicInterpolator reduce_y(const double y)

cdef class CySpline:
    cdef Spline *scpp

    def __init__(self, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] x, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] f):
        cdef vector[double] xvec = vector[double](len(x))
        cdef vector[double] fvec = vector[double](len(f))
        for i in range(len(x)):
            xvec[i] = x[i]
            fvec[i] = f[i]
        self.scpp = new Spline(xvec, fvec)

    def eval(self, double x):
        return self.scpp.evaluate(x)

    def deriv(self, double x):
        return self.scpp.derivative(x)

    def deriv2(self, double x):
        return self.scpp.derivative2(x)

cdef class CyCubicInterpolator:
    cdef CubicInterpolator *scpp

    def __init__(self, double x0, double dx, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] f):
        cdef vector[double] fvec = vector[double](len(f))
        for i in range(len(f)):
            fvec[i] = f[i]
        self.scpp = new CubicInterpolator(x0, dx, fvec)

    def eval(self, double x):
        return self.scpp.evaluate(x)

    def deriv(self, double x):
        return self.scpp.derivative(x)

    def deriv2(self, double x):
        return self.scpp.derivative2(x)

cdef class CySpline2D:
    cdef Spline2D *scpp

    def __init__(self, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] x, np.ndarray[ndim=1, dtype=np.float64_t, mode='c'] y, np.ndarray[ndim=2, dtype=np.float64_t, mode='c'] f):
        nx = x.shape[0]
        ny = y.shape[0]
        cdef vector[double] xvec = x
        cdef vector[double] yvec = y
        cdef vector[double] fvec = f.flatten()
        self.scpp = new Spline2D(xvec, yvec, fvec)

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
    

cdef class CyBicubicInterpolator:
    cdef BicubicInterpolator *scpp

    def __init__(self, double x0, double dx, int nx, double y0, double dy, int ny, np.ndarray[ndim=2, dtype=np.float64_t, mode='c'] f):
        cdef Matrix mz = Matrix(nx + 1, ny + 1)
        for i in range(nx + 1):
            for j in range(ny + 1):
                mz.set_value(i, j, f[i, j])
        self.scpp = new BicubicInterpolator(x0, dx, nx, y0, dy, ny, mz)

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

