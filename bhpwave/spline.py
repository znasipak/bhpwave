from bhpwaveformcy import CyCubicSpline, CyBicubicSpline, CyTricubicSpline
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
    A class for producing a cubic spline of a function :math:`f(x)` given its values
    :math:`f_{i} = f(x_i)` where :math:`x_i = x_0, x_1, \\dots , x_N` is a grid of :math:`(N+1)` uniformly-spaced
    points.
    
    Parameters
    ----------
    x : 1d-array[double]
        A uniformly-spaced grid of points
    f : 1d-array[double]
        Function values corresponding to the grid points x
    bc : str (optional)
        Boundary value method. Valid options include "natural", "not-a-knot", "clamped", and "E(3)"
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
        """
        The 2D array of spline coefficients with dimensions (nx, 4).
        Data are ordered so that the element at index (i, mx)
        returns the same value as coeffs(i, mx)

        Returns
        -------
        3d-array[double]
        """
        return np.array([[self.base.coefficient(i, j) for j in range(4)] for i in range(self.nx)])

    def coeff(self, i, mx):
        """
        Returns the spline coefficients :math:`c_{i}^{(m_x)}` defined by the
        spline :math:`f_{i}`

        .. math::
            f_{i}(x) = \\sum_{m_x=0}^3 c_{i}^{(m_x)}(x-x_i)^{m_x}

        Parameters
        ----------
        i : int
            Coefficient for the domain :math:`x_i \leq x \leq x_{i+1}`
        mx : int
            Coefficient weighting :math:`(x-x_i)^{m_x}` in the spline series
    
        Returns
        -------
        double
        """
        return self.base.coefficient(i, mx)

    def eval(self, x):
        """
        Evaluates the spline at the point x

        Parameters
        ----------
        x : double
            dependent parameter

        Returns
        -------
        double
        """
        if isinstance(x, np.ndarray):
            return np.array([self.base.eval(xi) for xi in x])
        return self.base.eval(x)
    
    def deriv(self, x):
        """
        Evaluates the derivative of the spline at the point x

        Parameters
        ----------
        x : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv(x)
    
    def deriv2(self, x):
        """
        Evaluates the second derivative of the spline at the point x

        Parameters
        ----------
        x : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv2(x)

    def __call__(self, x):
        """
        Evaluates the spline at the point x

        Parameters
        ----------
        x : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.eval(x)
    
class BicubicSpline:
    """
    A class for producing a bicubic spline of a function :math:`f(x, y)` given its values
    :math:`f_{ij} = f(x_i, y_j)` where :math:`x_i = x_0, x_1, \\dots , x_N` is a grid of :math:`(N+1)` uniformly-spaced
    points and :math:`y_j = y_0, y_1, \\dots , y_M` is a grid of :math:`(M+1)` uniformly-spaced
    points. The input :math:`f_{ij}` is therefore structured as a :math:`(N+1) \\times (M+1)` matrix of function values 
    
    .. math::
        \\begin{align*}
        f(x_i, y_j) &= 
            \\begin{pmatrix}
                f_{00} & f_{01} & \\cdots & f_{0M}
                \\\\
                f_{10} & f_{11} & \\cdots & f_{1M}
                \\\\
                \\vdots  &  \\vdots & \\ddots & \\vdots
                \\\\
                f_{N0} & f_{N1} & \\cdots & f_{NM}
            \\end{pmatrix}
        \\end{align*}
    
    
    Parameters
    ----------
    x : 1d-array[double]
        A uniformly-spaced grid of points
    y : 1d-array[double]
        A uniformly-spaced grid of points
    f : 2d-array[double]
        Function values corresponding to the grid points x, y
    bc : str (optional)
        Boundary value method. Valid options include "natural", "not-a-knot", "clamped", and "E(3)"
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
        """
        Evaluates the spline at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.eval(x, y)

    def deriv_x(self, x, y):
        """
        Evaluates the partial derivative of the spline with respect to x at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_x(x, y)
    
    def deriv_y(self, x, y):
        """
        Evaluates the partial derivative of the spline with respect to y at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_y(x, y)
    
    def deriv_xx(self, x, y):
        """
        Evaluates the second partial derivative of the spline with respect to x at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_xx(x, y)
    
    def deriv_yy(self, x, y):
        """
        Evaluates the second partial derivative of the spline with respect to y at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_yy(x, y)
    
    def deriv_xy(self, x, y):
        """
        Evaluates the mixed partial derivative of the spline with respect to x and y at the point (x, y)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_xy(x, y)
    
    def coeff(self, i, j, mx, my):
        """
        Returns the spline coefficients :math:`c_{ij}^{(m_x,m_y)}` defined by the
        spline :math:`f_{ij}`

        ..math::
            f_{ij}(x, y) = \\sum_{m_x=0}^3 \\sum_{m_y=0}^3 c_{ij}^{(m_x,m_y)}(x-x_i)^{m_x}(y-y_j)^{m_y}

        Parameters
        ----------
        i : int
            Coefficient for the domain :math:`x_i \leq x \leq x_{i+1}`
        j : int
            Coefficient for the domain :math:`y_j \leq y \leq y_{j+1}`
        mx : int
            Coefficient weighting :math:`(x-x_i)^{m_x}` in the spline series
        my : int
            Coefficient weighting :math:`(y-y_j)^{m_y}` in the spline series
    
        Returns
        -------
        double
        """
        return self.base.coefficient(i, j, mx, my)

    def __call__(self, x, y):
        return self.base.eval(x, y)
    
class TricubicSpline:
    """
    A class for producing a bicubic spline of a function :math:`f(x, y, z)` given its values
    :math:`f_{ijk} = f(x_i, y_j, z_k)` where :math:`x_i = x_0, x_1, \\dots , x_N` is a grid of :math:`(N+1)` uniformly-spaced
    points, :math:`y_j = y_0, y_1, \\dots , y_M` is a grid of :math:`(M+1)` uniformly-spaced
    points, and :math:`z_k = z_0, z_1, \\dots , z_O` is a grid of :math:`(O+1)` uniformly-spaced
    points. The input :math:`f_{ijk}` is therefore structured as a :math:`(N+1) \\times (M+1) \\times (O+1)` tensor of function values 
    
    .. math::
        \\begin{align*}
        f(x_i, y_j,z_k) &= 
            \\begin{pmatrix}
                f_{00k} & f_{01k} & \\cdots & f_{0Mk}
                \\\\
                f_{10k} & f_{11k} & \\cdots & f_{1Mk}
                \\\\
                \\vdots  &  \\vdots & \\ddots & \\vdots
                \\\\
                f_{N0k} & f_{N1k} & \\cdots & f_{NMk}
            \\end{pmatrix}
        \\end{align*}

    where each entry is a vector of length :math:`O+1` in the :math:`z`-dimension
    
    Parameters
    ----------
    x : 1d-array[double]
        A uniformly-spaced grid of points
    y : 1d-array[double]
        A uniformly-spaced grid of points
    z : 1d-array[double]
        A uniformly-spaced grid of points
    f : 3d-array[double]
        Function values corresponding to the grid points x, y, z or pre-computed spline coefficients
    bc : str (optional)
        Boundary value method. Valid options include "natural", "not-a-knot", "clamped", and "E(3)"
    """
    def __init__(self, x, y, z, f, bc = "E(3)"):
        self.boundary_conditions_dict = cubic_spline_bc_dict
        self.available_boundary_conditions = self.boundary_conditions_dict.keys()
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(f, np.ndarray)
        assert ((x.shape[0], y.shape[0], z.shape[0]) == (f.shape[0], f.shape[1], f.shape[2]) or (x.shape[0] - 1, y.shape[0] - 1, 64*(z.shape[0] - 1)) == (f.shape[0], f.shape[1], f.shape[2])), "Shapes of arrays {}, {}, {}, and {} do not match".format(x.shape, y.shape, z.shape, f.shape)

        self.x0 = x[0]
        self.y0 = y[0]
        self.z0 = z[0]
        self.dx = x[1]-self.x0
        self.dy = y[1]-self.y0
        self.dz = z[1]-self.z0
        self.nx = x.shape[0] - 1
        self.ny = y.shape[0] - 1
        self.nz = z.shape[0] - 1

        dx_array = np.diff(x)
        dy_array = np.diff(y)
        dz_array = np.diff(z)

        assert np.allclose(dx_array, self.dx*np.ones(dx_array.shape[0])), "Sampling points in x are not evenly spaced"
        assert np.allclose(dy_array, self.dy*np.ones(dy_array.shape[0])), "Sampling points in y are not evenly spaced"
        assert np.allclose(dz_array, self.dz*np.ones(dz_array.shape[0])), "Sampling points in z are not evenly spaced"

        self.base = CyTricubicSpline(self.x0, self.dx, self.nx, self.y0, self.dy, self.ny, self.z0, self.dz, self.nz, np.ascontiguousarray(f), self.boundary_conditions_dict[bc])

    def check_boundary_conditions(self, method):
        if method not in self.available_boundary_conditions:
            raise ValueError("No available method " + method)

    def eval(self, x, y, z):
        """
        Evaluates the spline at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.eval(x, y, z)

    def deriv_x(self, x, y, z):
        """
        Evaluates the partial derivative of the spline with respect to x at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_x(x, y, z)
    
    def deriv_y(self, x, y, z):
        """
        Evaluates the partial derivative of the spline with respect to y at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_y(x, y, z)
    
    def deriv_z(self, x, y, z):
        """
        Evaluates the partial derivative of the spline with respect to z at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_z(x, y, z)
    
    def deriv_xx(self, x, y, z):
        """
        Evaluates the second partial derivative of the spline with respect to x at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_xx(x, y, z)
    
    def deriv_yy(self, x, y, z):
        """
        Evaluates the second partial derivative of the spline with respect to y at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_yy(x, y, z)
    
    def deriv_zz(self, x, y, z):
        """
        Evaluates the second partial derivative of the spline with respect to z at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_zz(x, y, z)
    
    def deriv_xy(self, x, y, z):
        """
        Evaluates the mixed partial derivative of the spline with respect to x and y at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_xy(x, y, z)
    
    def deriv_xz(self, x, y, z):
        """
        Evaluates the mixed partial derivative of the spline with respect to x and z at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_xz(x, y, z)
    
    def deriv_yz(self, x, y, z):
        """
        Evaluates the mixed partial derivative of the spline with respect to y and z at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.deriv_yz(x, y, z)
    
    def coeff(self, i, j, k, mx, my, mz):
        """
        Returns the spline coefficients :math:`c_{ijk}^{(m_x,m_y,m_z)}` defined by the
        spline :math:`f_{ijk}:math:`

        .. math::
            f_{ijk}(x, y, z) = \\sum_{m_x=0}^3 \\sum_{m_y=0}^3 \\sum_{m_z=0}^3 c_{ijk}^{(m_x,m_y,m_z)}(x-x_i)^{m_x}(y-y_j)^{m_y}(z-z_k)^(m_z)

        Parameters
        ----------
        i : int
            Coefficient for the domain :math:`x_i \leq x \leq x_{i+1}`
        j : int
            Coefficient for the domain :math:`y_j \leq y \leq y_{j+1}`
        k : int
            Coefficient for the domain :math:`z_k \leq z \leq z_{k+1}`
        mx : int
            Coefficient weighting :math:`(x-x_i)^{m_x}` in the spline series
        my : int
            Coefficient weighting :math:`(y-y_j)^{m_y}` in the spline series
        mz : int
            Coefficient weighting :math:`(z-z_k)^{m_z}` in the spline series
    
        Returns
        -------
        double
        """
        return self.base.coefficient(i, j, k, mx, my, mz)
    
    def coefficients(self):
        """
        The 3D array of spline coefficients with dimensions (nx, ny, 64*nz).
        Data are ordered so that the element at index (i, j, k, 4*(4*(4*k + mx) + my) + mz)
        returns the same value as coeffs(i, j, k, mx, my, mz)

        Returns
        -------
        3d-array[double]
        """
        return self.base.coefficients()

    def __call__(self, x, y, z):
        """
        Evaluates the spline at the point (x, y, z)

        Parameters
        ----------
        x : double
            dependent parameter
        y : double
            dependent parameter
        z : double
            dependent parameter

        Returns
        -------
        double
        """
        return self.base.eval(x, y, z)
