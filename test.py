#test.py

from bhpwave.spline import BicubicSpline, TricubicSpline
from bhpwaveformcy import CyTricubicSpline
import numpy as np

x = np.linspace(0.5, 1, 200)
y = np.linspace(-1, 1, 400)
z = np.linspace(2, 3, 500)
XYZ = np.meshgrid(y, x, z)

def test_func(x, y, z):
    return np.sin(x)*np.cos(y)/(1-z**2)

fvec = test_func(XYZ[1], XYZ[0], XYZ[2])

spl = TricubicSpline(x, y, z, fvec)

print(spl(0.5, -1, 2) , test_func(0.5, -1, 2), fvec[0, 0, 0])