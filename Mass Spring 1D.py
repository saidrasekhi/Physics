'''
Program to solve 1D mass and spring system in presence of friction

Equations

F = dp/dt = -kx

'''

import pylab
import numpy
import scipy

from pylab import *
from numpy import *
from scipy import *

#
# initial condition, constants
#
h0 = .50     # initial deviation from equilibrium in m
v0 = 1.0     # m/s
t0 = 0.0     # s

#
# constants
#
time = 10.        # duration in sec
dt = 1.0e-3       # time step in sec
nstep = int(time/dt) # number of time steps

k = 10.0         # spring constant in N/m
m = 1.0          # mass in kg
k1_drag = 0.5    # kg/s

eta = 0.8        # v_in/v_out in an inellastic encounter
Energy = 0.5 * m * pow(v0,2) + 0.5 * k * pow(h0,2)

#
# array definitions
#
h = zeros((nstep,), dtype=float32)  # height in meters
v = zeros((nstep,), dtype=float32)  # velocity in m/sec
a = zeros((nstep,), dtype=float32)  # acceleration in m/sec^2
t = zeros((nstep,), dtype=float32)  # time in sec
K = zeros((nstep,), dtype=float32)  # kinetic energy
V = zeros((nstep,), dtype=float32)  # potential energy
En = zeros((nstep,), dtype=float32) # total energy

# initial conditions
i = 0
v[0] = v0
h[0] = h0
K[0] = 0.5 * m * pow(v0,2)
V[0] = 0.5 * k * pow(h0,2)
t[0] = t0
En[0] = K[0] + V[0]

for i in arange(1,nstep):   # main loop

    a[i-1] = -(k/m) * h[i-1] - (k1_drag/m) * v[i-1]
    v[i] = v[i-1] + a[i-1] * dt
    dE = - k1_drag * abs(v[i-1]) * abs(v[i-1]) * dt # drag
    Energy = Energy + dE
    h[i] = h[i-1] + v[i-1]*dt
    t[i] = t[i-1] + dt
    K[i] = 0.5 * m * pow(v[i],2)
    V[i] = 0.5 * k * pow(h[i],2)
    En[i] = K[i] + V[i]

tz = [0, time]
hz = [0, 0]
#
pylab.figure()
plot(t, h, tz, hz,'k')
xlabel('Time (s)', fontsize=18)
ylabel('Amplitude (m)', fontsize=18)
tex = r'$F\,=\,-kx\, -\, k_dv$'
text(5, .4, tex, fontsize=24)
savefig('spring_oscil.png', dpi=100)
#
pylab.figure()
plot(t, K, t, V, t, En, 'r', tz, hz,'k')
xlabel('Time (s)', fontsize=18)
ylabel('Energy (J)', fontsize=18)
legend(('Kinetic', 'Potential', 'Total'),'upper right', shadow=True)
savefig('spring_energy.png', dpi=100)
show()
