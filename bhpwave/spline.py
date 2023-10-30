from bhpwaveformcy import CyCubicSpline, CyBicubicSpline
import numpy as np

cubic_spline_bc_dict = {
    "natural": 0,
    "not-a-knot": 1,
    "clamped": 2,
    "E(3)": 3,
    "natural-alt": 4 
}

class CubicSpline:
    """
    A class for producing a cubic spline of a function f(x) given a grid of N+1 uniformly-spaced
    points x_i = x_0, x_1, ... , x_N and function values f(x_i) = f_0, f_1, ..., f_N
    
    :param x: A uniformly-spaced grid of points
    :type x: 1d-array[double]
    :param f: Function values corresponding to the grid points x
    :type f: 1d-array[double]
    :param method: Boundary value method
    :type f: str
    """
    def __init__(self, x, f, bc = "E(3)"):
        self.boundary_conditions_dict = cubic_spline_bc_dict
        self.available_boundary_conditions = self.boundary_conditions_dict.keys()

        assert isinstance(x, np.ndarray)
        assert isinstance(f, np.ndarray)
        assert x.shape == f.shape, "Shapes of arrays {} and {} do not match".format(x.shape, f.shape)

        self.x0 = x[0]
        self.dx = x[1] - self.x0
        self.nx = f.shape[0] - 1

        dx_array = x[1:] - x[:-1]
        assert np.allclose(dx_array, self.dx*np.ones(dx_array.shape[0])), "Sampling points are not evenly spaced"
        self.check_boundary_conditions(bc)
        
        self.base = CyCubicSpline(self.x0, self.dx, np.ascontiguousarray(f), self.boundary_conditions_dict[bc])
    
    def check_boundary_conditions(self, method):
        if method not in self.available_boundary_conditions:
            raise ValueError("No available method " + method)

    @property
    def coefficients(self):
        return np.array([[self.base.coefficient(i, j) for j in range(4)] for i in range(self.nx)])

    def coefficient(self, i, j):
        return self.base.coefficient(i, j)

    def eval(self, x):
        if isinstance(x, np.ndarray):
            return np.array([self.base.eval(xi) for xi in x])
        return self.base.eval(x)
    
    def deriv(self, x):
        return self.base.deriv(x)
    
    def deriv2(self, x):
        return self.base.deriv2(x)

    def __call__(self, x):
        return self.eval(x)
    
class BicubicSpline:
    """
    A class for producing a bicubic spline of a function f(x, y) given a grid of (N+1) uniformly-spaced
    points x_i = x_0, x_1, ... , x_N, a grid of (M+1) uniformly-spaced
    points y_j = y_0, y_1, ... , y_M and (N+1) x (M+1) matrix of function values 
    f(x_i, y_j) = f_{00}, f_{01}, ...  , f_{0M},
                  f_{10}, f_{11}, ...  , f_{1M},
                   ...  ,  ...  , ...  , ...  ,
                  f_{N0}, f_{N1}, ...  , f_{NM},
    
    :param x: A uniformly-spaced grid of points
    :type x: 1d-array[double]
    :param y: A uniformly-spaced grid of points
    :type y: 1d-array[double]
    :param f: Function values corresponding to the grid points x, y
    :type f: 2d-array[double]
    """
    def __init__(self, x, y, f, bc = "E(3)"):
        self.boundary_conditions_dict = cubic_spline_bc_dict
        self.available_boundary_conditions = self.boundary_conditions_dict.keys()
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(f, np.ndarray)
        assert (x.shape[0], y.shape[0]) == (f.shape[0], f.shape[1]), "Shapes of arrays {}, {}, and {} do not match".format(x.shape, y.shape, f.shape)

        self.x0 = x[0]
        self.y0 = y[0]
        self.dx = x[1]-self.x0
        self.dy = y[1]-self.y0
        self.nx = f.shape[0] - 1
        self.ny = f.shape[1] - 1

        dx_array = x[1:] - x[:-1]
        dy_array = y[1:] - y[:-1]
        assert np.allclose(dx_array, self.dx*np.ones(dx_array.shape[0])), "Sampling points in x are not evenly spaced"
        assert np.allclose(dy_array, self.dy*np.ones(dy_array.shape[0])), "Sampling points in y are not evenly spaced"

        self.base = CyBicubicSpline(self.x0, self.dx, self.nx, self.y0, self.dy, self.ny, np.ascontiguousarray(f), self.boundary_conditions_dict[bc])

    def check_boundary_conditions(self, method):
        if method not in self.available_boundary_conditions:
            raise ValueError("No available method " + method)

    def eval(self, x, y):
        return self.base.eval(x, y)

    def deriv_x(self, x, y):
        return self.base.deriv_x(x, y)
    
    def deriv_y(self, x, y):
        return self.base.deriv_y(x, y)
    
    def deriv_xx(self, x, y):
        return self.base.deriv_xx(x, y)
    
    def deriv_yy(self, x, y):
        return self.base.deriv_yy(x, y)
    
    def deriv_xy(self, x, y):
        return self.base.deriv_xy(x, y)

    def __call__(self, x, y):
        return self.base.eval(x, y)
