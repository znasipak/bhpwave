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
        pass

    cdef cppclass CubicInterpolator:
        CubicInterpolator(double x0, double dx, const vector[double] &y)
        CubicInterpolator(const vector[double] &x, const vector[double] &y)

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

