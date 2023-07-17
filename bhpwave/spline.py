from bhpwaveformcy import CySpline, CyCubicInterpolator
from bhpwaveformcy import CySpline2D, CyBicubicInterpolator

class Spline:
    def __init__(self, x, f):
        self.base = CySpline(x, f)
    
    def deriv(self, x):
        return self.base.deriv(x)
    
    def deriv2(self, x):
        return self.base.deriv2(x)

    def __call__(self, x):
        return self.base.eval(x)
    
class CubicInterpolant:
    def __init__(self, x0, dx, f):
        self.base = CyCubicInterpolator(x0, dx, f)

    def deriv(self, x):
        return self.base.deriv(x)
    
    def deriv2(self, x):
        return self.base.deriv2(x)

    def __call__(self, x):
        return self.base.eval(x)
    

class Spline2D:
    def __init__(self, x, y, f):
        self.base = CySpline2D(x, y, f)

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
    
class BicubicInterpolant:
    def __init__(self, x, y, f):
        x0 = x[0]
        y0 = y[0]
        dx = x[1]-x0
        dy = y[1]-y0
        nx = f.shape[0] - 1
        ny = f.shape[1] - 1
        self.base = CyBicubicInterpolator(x0, dx, nx, y0, dy, ny, f)

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
