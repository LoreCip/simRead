import numpy as np
from numba import njit

from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
import scipy.signal as sg

def find_zeros(x, arr):

    out = []
    
    interp = CubicSpline(x, arr)
    pmask = np.sign(arr)
    old = pmask[0]
    
    for i, sign in enumerate(pmask[1:]):
        if sign == old: continue
        root = brentq(interp, x[i], x[i+1])
        out.append([root, interp(root)])
        old = sign
    
    return np.array(out)

@njit
def expan(x, B, E):

    theta = 1 - x**2 * np.sin(2 * B / x**2)**2 / ( 4 * (1 + E) )
    im = np.argmin(theta)
    theta[im] = theta[im + 1] 
    
    return theta

@njit
def fxb(x, alpha, rs):
    return x**3 + alpha * x - rs


def find_i(idx, time, sim):
    # Check distance from wanted time
    dt = np.abs(time - sim.get(idx, 't')) 

    # Going backwards
    for i in range(idx-1, 0, -1):
        # Find new time
        t = sim.get(i, 't')
        # Check its distance from wanted time
        dt_tmp = np.abs(time - t)
        # If it is getting closer save it
        if dt_tmp < dt:
            dt = dt_tmp
        # If it is going further stop the cycle
        elif dt_tmp > dt:
            i = i + 1
            break
    return i
