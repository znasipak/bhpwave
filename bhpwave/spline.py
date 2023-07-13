from bhpwaveformcy import CySpline, CyCubicInterpolator

class Spline:
    def __init__(self, x, f):
        self.base = CySpline(x, f)

    def __call__(self, x):
        return self.base.eval(x)
    
class CubicInterpolant:
    def __init__(self, x0, dx, f):
        self.base = CyCubicInterpolator(x0, dx, f)

    def __call__(self, x):
        return self.base.eval(x)
