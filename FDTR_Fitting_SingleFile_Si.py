import math
import cmath
import tkinter
from tkinter import Tk  # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename
from scipy.optimize import leastsq
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import numpy
import numpy as np

# This program can be used when you are using two photodetectors simultaneously and all phase lag data is contained
# in one file.

deleteEnd = 1  # Deletes defined number of initial points from data
deleteBeg = 1  # Deletes defined number of end points from data

w0 = 6.9 * (10 ** (-6))  # Pump Diameter in m
w1 = 6.4 * (10 ** (-6))  # Probe Diameter in m
w0 /= 2  # Converts diameter to radius
w1 /= 2  # Converts diameter to radius

w = np.array([w0, w1], dtype=float)  # Creates an array with both pump and probe diameters

x0 = np.array([145, 1 * (10 ** 7)], dtype=float)  # Initial guesses for unknown parameters

# Fitting function. Note that "cos" is NOT cosine. It's the unknown parameter(s). Python is 0 indexed, so cos[0] is
# the first unknown parameter


def fun(frequency, cos, wa):
    w0a = wa[0]
    w1a = wa[1]
    kr = [220, 0, cos[0]]  # Vector of in-plane thermal conductivities for each layer; even layers are 0 always
    # because they are interfaces.
    kz = [220, cos[1], cos[0]]  # Vector of through-plane thermal conductivities (odd layers) and thermal boundary
    # conductances (even layers)
    cv = [2.48 * (10 ** 6), 0, 1.8 * (10 ** 6)]  # Vector of volumetric heat capacities for each layer, which is the product of specific heat capacity (cp) and density (rho). Even layers are always 0 because they are interfaces.
    t = [143 * (10 ** (-9)), 0, 1]  # Thickness of each layer; even layers are always 0 because they are interfaces.
    k = numpy.logspace(2.0, 9.0, num=100)  # Integration parameter for Hankel Transform
    omega = np.array(frequency * math.pi * 2)  # Defined frequency vector in radians
    phi = np.ones((len(omega),), dtype=float)  # Defining the phase vector for faster evaluation in for loop.
    phi_f = np.ones((len(omega),), dtype=float)
    for f in range(len(omega)):
        h = hankel(k, wa, omega[f], kr, kz, cv, t)
        hint = (1 / (2 * math.pi)) * numpy.trapz(h, x=k)
        phi_f = math.atan(hint.imag / hint.real)
        phi[f] = (180/math.pi) * phi_f
    return phi


def hankel(k, wa, omega, kr, kz, cv, t):  #Hankel transform for infinite mulilayer solution; developed by R. Warzoha
    w0a = wa[0]
    w1a = wa[1]
    hint = np.ones((len(k),), dtype='complex_')  # np.ones((len(k),), dtype=float)
    q = np.ones((len(t),), dtype='complex_')
    nn = [[], [], []]  #np.array([[[1, 1], [1, 1]], [[1, 1], [1, 1]], [[1, 1], [1, 1]]], dtype='complex_')
    for i in range(len(k)):
        mm = [[1, 0],
              [0, 1]]
        for p in range(len(t)-1, -1, -1):
            q[p] = ((kr[p] * (k[i] ** 2) + cv[p] * 1j * omega) / (kz[p])) ** 0.5
            if t[p] != 0:
                nn = [[1, -(cmath.tanh(q[p] * t[p])) / (kz[p] * q[p])], [-kz[p] * q[p] * cmath.tanh(q[p] * t[p]), 1]]
                mm = np.dot(mm, nn)
            else:
                nn = [[1, -kz[p] ** (-1)], [0, 1]]
                mm = np.dot(mm, nn)
        ndc = -mm[1, 1] / mm[1, 0]
        hint[i] = k[i] * ndc * cmath.exp((-(k[i] ** 2) * (w0a ** 2 + w1a ** 2)) / 8)
    return hint


def residuals(cos, y, frequency, wa):
    return y - fun(frequency, cos, wa)  # Minimization of SSE


root = tkinter.Tk()
Tk().withdraw()  # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename()  # show an "Open" dialog box and return the path to the selected file
root.destroy()

data = np.loadtxt(filename, delimiter=';', skiprows=6, dtype=float)  # Loads the text file from Zurich; this is the signal file


freq = data[deleteBeg:-deleteEnd, 0]  # Pulls the measured frequency data
phi_data = data[deleteBeg:-deleteEnd, 1]  # Pulls the measured phase data
phi_ref = data[deleteBeg:-deleteEnd, 2]  # Pulls the reference phase data
phi_exp = phi_data - phi_ref  # Calculates the phase lag between pump and probe signals

# Phase lag correction if lock-in is recording 180 degrees out of phase, which happens randomly:
if phi_exp.any() < -180:
    phi_exp += 180
elif phi_exp.any() > 0:
    phi_exp -= 180


# Least squares fitting routine:
gg = least_squares(residuals, x0, loss='soft_l1', args=(phi_exp, freq, w), bounds=(0.00005, np.inf))
print(gg.x) # Print the results of the fitting


# Plot the results:
fig, ax = plt.subplots(figsize=(5, 3))
ax.scatter(x=freq, y=phi_exp, marker='o', c='r', edgecolor='b')
ax.plot(freq, fun(freq, gg.x, w))
ax.set_xscale('log')
ax.set_title('Raw Data')
ax.set_ylabel('Phase Lag [deg]')
ax.set_xlabel('Modulation Frequency [Hz]')
ax.set_xlim(xmin=freq[0], xmax=freq[-1] + 100)
fig.tight_layout()
plt.show()
