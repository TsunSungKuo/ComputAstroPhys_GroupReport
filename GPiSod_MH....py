import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import minimize

#--------------------------------------------------------------------
# parameters
#--------------------------------------------------------------------
# constants
L        = 1.0       # 1-D computational domain size
N_In     = 128       # number of computing cells
cfl      = 1.0       # Courant factor
nghost   = 2         # number of ghost zones
gamma    = 5.0/3.0   # ratio of specific heats
end_time = 0.005       # simulation time

# derived constants
N  = N_In + 2*nghost    # total number of cells including ghost zones
dx = L/N_In             # spatial resolution

# plotting parameters
nstep_per_image = 1     # plotting frequency


# -------------------------------------------------------------------
# define initial condition
# -------------------------------------------------------------------
def InitialCondition( x ):
    if ( x < 0.5*L ):
        d = 1.0  # density
        u = 0.0  # velocity x
        v = 0.0  # velocity y
        w = 0.0  # velocity z
        P = 1.0  # pressure
        E = P/(gamma-1.0) + 0.5*d*( u**2.0 + v**2.0 + w**2.0 )   # energy density
    else:
        d = 0.125
        u = 0.0
        v = 0.0
        w = 0.0
        P = 0.1
        E = P/(gamma-1.0) + 0.5*d*( u**2.0 + v**2.0 + w**2.0 )
    return np.array([d, d*u, d*v, d*w, E])

# -------------------------------------------------------------------
# define boundary condition by setting ghost zones
# -------------------------------------------------------------------
def BoundaryCondition(U):
    U[0:nghost]   = U[nghost]
    U[N-nghost:N] = U[N-nghost-1]

# -------------------------------------------------------------------
# compute pressure
# -------------------------------------------------------------------
def ComputePressure(d, px, py, pz, e):
    P = (gamma-1.0) * (e - 0.5*(px**2.0 + py**2.0 + pz**2.0)/d)
    epsilon = 1e-8
    P = np.maximum(P, epsilon)  # Ensure pressure is always positive
    return P

# -------------------------------------------------------------------
# compute time-step by the CFL condition
# -------------------------------------------------------------------
def ComputeTimestep(U):
    P = ComputePressure(U[:,0], U[:,1], U[:,2], U[:,3], U[:,4])
    a = (gamma * P / U[:,0])**0.5
    u = np.abs(U[:,1] / U[:,0])
    max_info_speed = np.amax(u + a)
    dt_cfl = cfl * dx / max_info_speed
    dt_end = end_time - t
    return min(dt_cfl, dt_end)

# -------------------------------------------------------------------
# Define the RBF (Squared Exponential) kernel function
# -------------------------------------------------------------------
def rbf_kernel(X1, X2, length_scale, variance):
    sqdist = np.sum(X1 ** 2, 1).reshape(-1, 1) + np.sum(X2 ** 2, 1) - 2 * np.dot(X1, X2.T)
    return variance * np.exp(-0.5 / length_scale ** 2 * sqdist)

# Define the log marginal likelihood function
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
    initial_theta = [0.1, 1.0, 1e-8]
    result = minimize(log_marginal_likelihood, initial_theta, args=(X_train, y_train), bounds=((1e-5, None), (1e-5, None), (1e-5, None)))
    optimized_theta = result.x
    length_scale, variance, noise_variance = optimized_theta
    
    K = rbf_kernel(X_train, X_train, length_scale, variance) + noise_variance * np.eye(len(X_train))
    L = np.linalg.cholesky(K)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))
    
    K_s = rbf_kernel(X_train, X_test, length_scale, variance)
    K_ss = rbf_kernel(X_test, X_test, length_scale, variance) + noise_variance * np.eye(len(X_test))
    
    mu_s = K_s.T.dot(alpha)
    v = np.linalg.solve(L, K_s)
    cov_s = K_ss - v.T.dot(v)
    sigma_s = np.sqrt(np.diag(cov_s))
    
    return mu_s, sigma_s, optimized_theta

# -------------------------------------------------------------------
# GP-based data reconstruction
# -------------------------------------------------------------------
def GPReconstruction(U):
    x_train = np.linspace(0, L, N_In).reshape(-1, 1)
    U_gp = np.empty(U.shape)
    for i in range(U.shape[1]):  # for each conserved variable
        y_train = U[nghost:N-nghost, i]
        mu, _, _ = train_gp_and_predict(x_train, y_train, x_train)
        U_gp[nghost:N-nghost, i] = mu
    return U_gp

# -------------------------------------------------------------------
# convert conserved variables to fluxes
# -------------------------------------------------------------------
def Conserved2Flux(U):
    flux = np.empty(5)
    P = ComputePressure(U[0], U[1], U[2], U[3], U[4])
    u = U[1] / U[0]
    flux[0] = U[1]
    flux[1] = u * U[1] + P
    flux[2] = u * U[2]
    flux[3] = u * U[3]
    flux[4] = u * (U[4] + P)
    return flux

# -------------------------------------------------------------------
# Roe's Riemann solver
# -------------------------------------------------------------------
def Roe(L, R):
    P_L = ComputePressure(L[0], L[1], L[2], L[3], L[4])
    P_R = ComputePressure(R[0], R[1], R[2], R[3], R[4])
    H_L = (L[4] + P_L) / L[0]
    H_R = (R[4] + P_R) / R[0]

    rhoL_sqrt = L[0]**0.5
    rhoR_sqrt = R[0]**0.5

    u  = (L[1]/rhoL_sqrt + R[1]/rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt)
    v  = (L[2]/rhoL_sqrt + R[2]/rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt)
    w  = (L[3]/rhoL_sqrt + R[3]/rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt)
    H  = (rhoL_sqrt*H_L  + rhoR_sqrt*H_R)  / (rhoL_sqrt + rhoR_sqrt)
    V2 = u*u + v*v + w*w

    a  = ((gamma-1.0)*(H - 0.5*V2))**0.5

    dU     = R - L
    amp    = np.empty(5)
    amp[2] = dU[2] - v*dU[0]
    amp[3] = dU[3] - w*dU[0]
    amp[1] = (gamma-1.0)/a**2.0 * (dU[0]*(H-u**2.0) + u*dU[1] - dU[4] + v*amp[2] + w*amp[3])
    amp[0] = 0.5/a * (dU[0]*(u+a) - dU[1] - a*amp[1])
    amp[4] = dU[0] - amp[0] - amp[1]

    EigenValue    = np.array([u-a, u, u, u, u+a])
    EigenVector_R = np.array([[1.0, u-a,   v,   w,  H-u*a],
                              [1.0,   u,   v,   w, 0.5*V2],
                              [0.0, 0.0, 1.0, 0.0,      v],
                              [0.0, 0.0, 0.0, 1.0,      w],
                              [1.0, u+a,   v,   w,  H+u*a]])

    flux_L = Conserved2Flux(L)
    flux_R = Conserved2Flux(R)

    amp *= np.abs(EigenValue)
    flux = 0.5*(flux_L + flux_R) - 0.5*amp.dot(EigenVector_R)
    return flux

# -------------------------------------------------------------------
# initialize animation
# -------------------------------------------------------------------
def init():
    line_d.set_xdata(x)
    line_u.set_xdata(x)
    line_p.set_xdata(x)
    return line_d, line_u, line_p

# -------------------------------------------------------------------
# update animation
# -------------------------------------------------------------------
def update(frame):
    global t, U

    if frame > 0:
        for step in range(nstep_per_image):
            BoundaryCondition(U)
            dt = ComputeTimestep(U)
            print("t = %13.7e --> %13.7e, dt = %13.7e" % (t, t+dt, dt))

            L, R = GPReconstruction(U), GPReconstruction(U)

            for j in range(1, N-1):
                flux_L = Conserved2Flux(L[j])
                flux_R = Conserved2Flux(R[j])
                dflux  = 0.5*dt/dx*(flux_R - flux_L)
                L[j]  -= dflux
                R[j]  -= dflux

            flux = np.empty((N, 5))
            for j in range(nghost, N-nghost+1):
                flux[j] = Roe(R[j-1], L[j])

            U[nghost:N-nghost] -= dt/dx*(flux[nghost+1:N-nghost+1] - flux[nghost:N-nghost])

            t = t + dt
            if (t >= end_time):
                anim.event_source.stop()
                break

    d = U[nghost:N-nghost,0]
    u = U[nghost:N-nghost,1] / U[nghost:N-nghost,0]
    P = ComputePressure(U[nghost:N-nghost,0], U[nghost:N-nghost,1], U[nghost:N-nghost,2], U[nghost:N-nghost,3], U[nghost:N-nghost,4])
    line_d.set_ydata(d)
    line_u.set_ydata(u)
    line_p.set_ydata(P)
    ax[0].set_title('t = %6.3f' % t)
    return line_d, line_u, line_p

# Set initial condition
t = 0.0
x = np.empty(N_In)
U = np.empty((N, 5))
for j in range(N_In):
    x[j] = (j + 0.5) * dx
    U[j + nghost] = InitialCondition(x[j])

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
ax[2].set_ylim( +0.0, 1.2 )

# Create movie
nframe = 99999999  # arbitrarily large
anim = animation.FuncAnimation(fig, func=update, init_func=init, frames=nframe, interval=200, repeat=False)
plt.show()
