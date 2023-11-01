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

    cdef cppclass CubicSpline:
        CubicSpline(double x0, double dx, const vector[double] &y, int method) except +
        CubicSpline(const vector[double] &x, const vector[double] &y, int method) except +

        double getSplineCoefficient(int i, int j)

        double evaluate(const double x)
        double derivative(const double x)
        double derivative2(const double x)

    cdef cppclass BicubicSpline:
        BicubicSpline(const vector[double] &x, const vector[double] &y, const Matrix &z, int method)
        BicubicSpline(double x0, double dx, int nx, double y0, double dy, int ny, const Matrix &z, int method)
        
        double evaluate(const double x, const double y)
        double derivative_x(const double x, const double y)
        double derivative_y(const double x, const double y)
        double derivative_xy(const double x, const double y)
        double derivative_xx(const double x, const double y)
        double derivative_yy(const double x, const double y)
        CubicSpline reduce_x(const double x)
        CubicSpline reduce_y(const double y)

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

