import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import minimize

# Parameters
L = 1.0       # 1-D computational domain size
N_In = 128    # number of computing cells
cfl = 1.0     # Courant factor
nghost = 1    # number of ghost zones
gamma = 5.0 / 3.0  # ratio of specific heats
end_time = 0.1     # simulation time

N = N_In + 2 * nghost  # total number of cells including ghost zones
dx = L / N_In          # spatial resolution

nstep_per_image = 1  # plotting frequency

# Initial condition
def InitialCondition(x):
    if x < 0.5 * L:
        d, u, P = 1.0, 0.0, 1.0
    else:
        d, u, P = 0.125, 0.0, 0.1
    E = P / (gamma - 1.0) + 0.5 * d * u**2.0
    return np.array([d, d*u, E])

# Boundary condition
def BoundaryCondition(U):
    U[0:nghost] = U[nghost]
    U[N-nghost:N] = U[N-nghost-1]

# Compute pressure
def ComputePressure(d, px, e):
    P = (gamma - 1.0) * (e - 0.5 * px**2.0 / d)
    return P

# Compute timestep
def ComputeTimestep(U):
    P = ComputePressure(U[:, 0], U[:, 1], U[:, 2])
    a = (gamma * P / U[:, 0])**0.5
    u = np.abs(U[:, 1] / U[:, 0])
    max_info_speed = np.amax(u + a)
    dt_cfl = cfl * dx / max_info_speed
    dt_end = end_time - t
    return min(dt_cfl, dt_end)

# Gaussian Process kernel function
def rbf_kernel(X1, X2, length_scale, variance):
    sqdist = np.sum(X1 ** 2, 1).reshape(-1, 1) + np.sum(X2 ** 2, 1) - 2 * np.dot(X1, X2.T)
    return variance * np.exp(-0.5 / length_scale ** 2 * sqdist)

# Log marginal likelihood function
def log_marginal_likelihood(theta, X, y):
    length_scale, variance, noise_variance = theta
    K = rbf_kernel(X, X, length_scale, variance) + noise_variance * np.eye(len(X))
    L = np.linalg.cholesky(K)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y))
    log_likelihood = -0.5 * np.dot(y.T, alpha)
    log_likelihood -= np.sum(np.log(np.diag(L)))
    log_likelihood -= len(X) / 2 * np.log(2 * np.pi)
    return -log_likelihood  # return the negative log likelihood for minimization

# Function to train GP and make predictions
def train_gp_and_predict(X_train, y_train, X_test):
    result = minimize(log_marginal_likelihood, initial_theta, args=(X_train, y_train), bounds=((1e-5, None), (1e-5, None), (1e-5, None)))
    optimized_theta = result.x
    length_scale, variance, noise_variance = optimized_theta
    
    K = rbf_kernel(X_train, X_train, length_scale, variance) + noise_variance * np.eye(len(X_train))
    L = np.linalg.cholesky(K)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))
    
    K_s = rbf_kernel(X_train, X_test, length_scale, variance)
    mu_s = K_s.T.dot(alpha)
    
    return mu_s, optimized_theta

# Convert conserved variables to fluxes
def Conserved2Flux(U):
    flux = np.empty(3)
    P = ComputePressure(U[0], U[1], U[2])
    u = U[1] / U[0]
    flux[0] = U[1]
    flux[1] = u * U[1] + P
    flux[2] = u * (U[2] + P)
    return flux

# Initialize animation
def init():
    line_d.set_xdata(x)
    line_u.set_xdata(x)
    line_p.set_xdata(x)
    return line_d, line_u, line_p

# Update animation
def update(frame):
    global t, d, u, P

    if frame > 0:
        for step in range(nstep_per_image):
            x_train = x.reshape(-1, 1)

            BoundaryCondition(U)

            # Estimate time-step from the CFL condition
            dt = ComputeTimestep(U)
            print("t = %13.7e --> %13.7e, dt = %13.7e" % (t, t + dt, dt))

            # Compute fluxes
            flux = np.empty((N, 3))
            for j in range(nghost, N-nghost+1):
                flux[j] = 0.5 * (Conserved2Flux(U[j]) + Conserved2Flux(U[j-1]) - dx/dt * (U[j] - U[j-1]))

            # Update the volume-averaged input variables by dt
            U[nghost:N-nghost] -= dt / dx * (flux[nghost+1:N-nghost+1] - flux[nghost:N-nghost])

            y_train_density = d
            y_train_velocity = u
            y_train_pressure = P
            
            d, _ = train_gp_and_predict(x_train, y_train_density, x_train)
            u, _ = train_gp_and_predict(x_train, y_train_velocity, x_train)
            P, _ = train_gp_and_predict(x_train, y_train_pressure, x_train)
            
            t += dt
            if t >= end_time:
                anim.event_source.stop()
                break

    line_d.set_ydata(d)
    line_u.set_ydata(u)
    line_p.set_ydata(P)
    ax[0].set_title('t = %6.3f' % t)
    return line_d, line_u, line_p

# Set initial condition
t = 0.0
x = np.linspace(0, L, N)
U = np.empty((N, 3))
for j in range(N):
    U[j] = InitialCondition((j + 0.5 - nghost) * dx)

d = U[:, 0]
u = U[:, 1]
P = U[:, 2]

# Initial hyperparameters for GP
initial_theta = [0.1, 1.0, 1e-8]

# Create figure
fig, ax = plt.subplots(3, 1, sharex=True, sharey=False, dpi=140)
fig.subplots_adjust(hspace=0.1, wspace=0.0)
line_d, = ax[0].plot([], [], 'r-o', ls='-', markeredgecolor='k', markersize=3)
line_u, = ax[1].plot([], [], 'b-o', ls='-', markeredgecolor='k', markersize=3)
line_p, = ax[2].plot([], [], 'g-o', ls='-', markeredgecolor='k', markersize=3)
ax[2].set_xlabel('x')
ax[0].set_ylabel('Density')
ax[1].set_ylabel('Velocity')
ax[2].set_ylabel('Pressure')
ax[0].set_xlim( 0.30, L-0.25 )
ax[0].set_ylim( +0.0, 1.2 )
ax[1].set_ylim( -0.2, 1.0 )
ax[2].set_ylim( +0.0, 1.8 )

# Create movie
nframe = 99999999
anim = animation.FuncAnimation(fig, func=update, init_func=init, frames=nframe, interval=10, repeat=False)
plt.show()
