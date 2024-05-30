import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def compute_smoothness_indicators(v):  #beta
    B = np.zeros(3)
    B[0] = (13/12)*(v[0] - 2*v[1] + v[2])**2 + (1/4)*(v[0] - 4*v[1] + 3*v[2])**2
    B[1] = (13/12)*(v[1] - 2*v[2] + v[3])**2 + (1/4)*(v[1] - v[3])**2
    B[2] = (13/12)*(v[2] - 2*v[3] + v[4])**2 + (1/4)*(3*v[2] - 4*v[3] + v[4])**2
    return B

def compute_weights(B, eps=1e-6):
    d = np.array([0.1, 0.6, 0.3])  # posituve linear weight??
    alpha = d / (eps + B)**2
    w = alpha / np.sum(alpha)
    return w

def weno_flux(v):
    f = np.zeros(len(v) - 4)
    for i in range(2, len(v) - 2):
        stencil = v[i-2:i+3]
        B = compute_smoothness_indicators(stencil)
        w = compute_weights(B)

        f_stencils = np.array([
            stencil[0]*1/3 - stencil[1]*7/6 + stencil[2]*11/6,
            -stencil[1]*1/6 + stencil[2]*5/6 + stencil[3]*1/3,
            stencil[2]*1/3 + stencil[3]*5/6 - stencil[4]*1/6
        ])

        f[i-2] = np.sum(w * f_stencils)
    return f

# Example usage:
u = np.sin(np.linspace(0, 2*np.pi, 100))  # Example initial condition
flux = weno_flux(u)

plt.figure(figsize=(10, 6))
plt.plot(u, label='Initial Condition (u)')
plt.plot(flux, label='Flux')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.title('WENO Flux Calculation')
plt.show()