import numpy as np
import pylab
kappa = 1
delta = 10.
a = 0.5
b = 1.0
k = np.pi / 2 / delta

def u(x, t):
    return  a / k * np.cos(k * x) * np.exp(-k * kappa * t) + a * (x - delta) + b

N = 1000
x = np.arange(N) / float(N - 1) * delta

pylab.plot(x, u(x, 0), x, u(x, 1), x, u(x, 5))

pylab.show()
