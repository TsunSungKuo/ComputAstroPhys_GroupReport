import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
L = 1.0  # 1-D computational domain size
N_In = 128  # number of computing cells
cfl = 1.0  # Courant factor
nghost = 2  # number of ghost zones
gamma = 5.0 / 3.0  # ratio of specific heats
end_time = 0.4  # simulation time

# Derived constants
N = N_In + 2 * nghost  # total number of cells including ghost zones
dx = L / N_In  # spatial resolution

# Define initial condition
def InitialCondition(x):
    if x < 0.5 * L:
        d = 1.25e3  # density
        u = 0.0  # velocity x
        v = 0.0  # velocity y
        w = 0.0  # velocity z
        P = 5.0e2  # pressure
        E = P / (gamma - 1.0) + 0.5 * d * (u**2.0 + v**2.0 + w**2.0)  # energy density
    else:
        d = 1.25e2
        u = 0.0
        v = 0.0
        w = 0.0
        P = 5.0e0
        E = P / (gamma - 1.0) + 0.5 * d * (u**2.0 + v**2.0 + w**2.0)

    return np.array([d, d * u, d * v, d * w, E])

# Boundary condition
def BoundaryCondition(U):
    U[0:nghost] = U[nghost]
    U[N - nghost:N] = U[N - nghost - 1]

# Compute pressure
def ComputePressure( d, px, py, pz, e ):
   P = (gamma-1.0)*( e - 0.5*(px**2.0 + py**2.0 + pz**2.0)/d )
   assert np.all( P > 0 ), "negative pressure !!"
   return P

# Compute time-step by the CFL condition
def ComputeTimestep(U):
    P = ComputePressure(U[:, 0], U[:, 1], U[:, 2], U[:, 3], U[:, 4])
    a = (gamma * P / U[:, 0])**0.5
    u = np.abs(U[:, 1] / U[:, 0])
    v = np.abs(U[:, 2] / U[:, 0])
    w = np.abs(U[:, 3] / U[:, 0])

    max_info_speed = np.amax(u + a)
    dt_cfl = cfl * dx / max_info_speed
    dt_end = end_time - t

    return min(dt_cfl, dt_end)

# Compute smoothness indicators
def compute_smoothness_indicators(v):
    B = np.zeros(3)
    B[0] = (13 / 12) * (v[0] - 2 * v[1] + v[2])**2 + (1 / 4) * (v[0] - 4 * v[1] + 3 * v[2])**2
    B[1] = (13 / 12) * (v[1] - 2 * v[2] + v[3])**2 + (1 / 4) * (v[1] - v[3])**2
    B[2] = (13 / 12) * (v[2] - 2 * v[3] + v[4])**2 + (1 / 4) * (3 * v[2] - 4 * v[3] + v[4])**2
    return B

# Compute weights
def compute_weights(B, eps=1e-6):
    d = np.array([0.1, 0.6, 0.3])
    alpha = d / (eps + B)**2
    w = alpha / np.sum(alpha)
    return w

# WENO flux calculation
def weno_flux(v):
    f = np.zeros(len(v) - 4)
    for i in range(2, len(v) - 2):
        stencil = v[i - 2:i + 3]
        B = compute_smoothness_indicators(stencil)
        w = compute_weights(B)

        f_stencils = np.array([
            stencil[0] * 1 / 3 - stencil[1] * 7 / 6 + stencil[2] * 11 / 6,
            -stencil[1] * 1 / 6 + stencil[2] * 5 / 6 + stencil[3] * 1 / 3,
            stencil[2] * 1 / 3 + stencil[3] * 5 / 6 - stencil[4] * 1 / 6
        ])

        f[i - 2] = np.sum(w * f_stencils)
    return f

# Data reconstruction using WENO
def DataReconstruction_WENO(U):
    W = np.empty((N, 5))
    L = np.empty((N, 5))
    R = np.empty((N, 5))

    for j in range(N):
        W[j] = Conserved2Primitive(U[j])

    for j in range(2, N - 2):
        for k in range(5):
            stencils = W[j - 2:j + 3, k]
            f = weno_flux(stencils)
            L[j, k] = f[0]
            R[j, k] = f[-1]

        L[j] = Primitive2Conserved(L[j])
        R[j] = Primitive2Conserved(R[j])

    return L, R

# Convert conserved variables to primitive variables
def Conserved2Primitive(U):
    W = np.empty(5)
    W[0] = U[0]
    W[1] = U[1] / U[0] if U[0] != 0 else 0
    W[2] = U[2] / U[0] if U[0] != 0 else 0
    W[3] = U[3] / U[0] if U[0] != 0 else 0
    W[4] = ComputePressure(U[0], U[1], U[2], U[3], U[4])
    return W

# Convert primitive variables to conserved variables
def Primitive2Conserved(W):
    U = np.empty(5)
    U[0] = W[0]
    U[1] = W[0] * W[1]
    U[2] = W[0] * W[2]
    U[3] = W[0] * W[3]
    U[4] = W[4] / (gamma - 1.0) + 0.5 * W[0] * (W[1]**2.0 + W[2]**2.0 + W[3]**2.0)
    return U

# Convert conserved variables to fluxes
def Conserved2Flux(U):
    flux = np.empty(5)
    P = ComputePressure(U[0], U[1], U[2], U[3], U[4])
    u = U[1] / U[0] if U[0] != 0 else 0
    flux[0] = U[1]
    flux[1] = u * U[1] + P
    flux[2] = u * U[2]
    flux[3] = u * U[3]
    flux[4] = u * (U[4] + P)
    return flux

# Roe's Riemann solver
def Roe( L, R ):
#  compute the enthalpy of the left and right states: H = (E+P)/rho
   P_L = ComputePressure( L[0], L[1], L[2], L[3], L[4] )
   P_R = ComputePressure( R[0], R[1], R[2], R[3], R[4] )
   H_L = ( L[4] + P_L )/L[0]
   H_R = ( R[4] + P_R )/R[0]

#  compute Roe average values
   rhoL_sqrt = L[0]**0.5
   rhoR_sqrt = R[0]**0.5

   u  = ( L[1]/rhoL_sqrt + R[1]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt )
   v  = ( L[2]/rhoL_sqrt + R[2]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt )
   w  = ( L[3]/rhoL_sqrt + R[3]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt )
   H  = ( rhoL_sqrt*H_L  + rhoR_sqrt*H_R  ) / ( rhoL_sqrt + rhoR_sqrt )
   V2 = u*u + v*v + w*w
#  check negative pressure
   assert H-0.5*V2 > 0.0, "negative pressure!"
   a  = ( (gamma-1.0)*(H - 0.5*V2) )**0.5

#  compute the amplitudes of different characteristic waves
   dU     = R - L
   amp    = np.empty( 5 )
   amp[2] = dU[2] - v*dU[0]
   amp[3] = dU[3] - w*dU[0]
   amp[1] = (gamma-1.0)/a**2.0 \
            *( dU[0]*(H-u**2.0) + u*dU[1] - dU[4] + v*amp[2] + w*amp[3] )
   amp[0] = 0.5/a*( dU[0]*(u+a) - dU[1] - a*amp[1] )
   amp[4] = dU[0] - amp[0] - amp[1]

#  compute the eigenvalues and right eigenvector matrix
   EigenValue    = np.array( [u-a, u, u, u, u+a] )
   EigenVector_R = np.array(  [ [1.0, u-a,   v,   w,  H-u*a],
                                [1.0,   u,   v,   w, 0.5*V2],
                                [0.0, 0.0, 1.0, 0.0,      v],
                                [0.0, 0.0, 0.0, 1.0,      w],
                                [1.0, u+a,   v,   w,  H+u*a] ]  )

#  compute the fluxes of the left and right states
   flux_L = Conserved2Flux( L )
   flux_R = Conserved2Flux( R )

#  compute the Roe flux
   amp *= np.abs( EigenValue )
   flux = 0.5*( flux_L + flux_R ) - 0.5*amp.dot( EigenVector_R )

   return flux

# Initialization function for animation
def init():
    line_d.set_data([], [])
    line_u.set_data([], [])
    line_p.set_data([], [])
    return line_d, line_u, line_p

# Update function for animation
def update( frame ):
   global t, U

#  for frame==0, just plot the initial condition
   if frame > 0:
      for step in range( nstep_per_image ):

#        set the boundary conditions
         BoundaryCondition( U )

#        estimate time-step from the CFL condition
         dt = ComputeTimestep( U )
         print( "t = %13.7e --> %13.7e, dt = %13.7e" % (t,t+dt,dt) )

#        data reconstruction
         L, R = DataReconstruction_WENO( U )

#        update the face-centered variables by 0.5*dt
         for j in range( 1, N-1 ):
            flux_L = Conserved2Flux( L[j] )
            flux_R = Conserved2Flux( R[j] )
            dflux  = 0.5*dt/dx*( flux_R - flux_L )
            L[j]  -= dflux
            R[j]  -= dflux

#        compute fluxes
         flux = np.empty( (N,5) )
         for j in range( nghost, N-nghost+1 ):
#           R[j-1] is the LEFT state at the j+1/2 inteface
            flux[j] = Roe( R[j-1], L[j]  )

#        update the volume-averaged input variables by dt
         U[nghost:N-nghost] -= dt/dx*( flux[nghost+1:N-nghost+1] - flux[nghost:N-nghost] )

#        update time
         t = t + dt
         if ( t >= end_time ):
            anim.event_source.stop()
            break

#  plot
   d = U[nghost:N-nghost,0]
   u = U[nghost:N-nghost,1] / U[nghost:N-nghost,0]
   P = ComputePressure( U[nghost:N-nghost,0], U[nghost:N-nghost,1], U[nghost:N-nghost,2],
                        U[nghost:N-nghost,3], U[nghost:N-nghost,4] )
   line_d.set_ydata( d )
   line_u.set_ydata( u )
   line_p.set_ydata( P )
#  ax[0].legend( loc='upper right', fontsize=12 )
#  ax[1].legend( loc='upper right', fontsize=12 )
#  ax[2].legend( loc='upper right', fontsize=12 )
   ax[0].set_title( 't = %6.3f' % (t) )

   return line_d, line_u, line_p



# Main simulation loop
t = 0.0
U = np.empty((N, 5))

x = np.linspace(-nghost * dx + 0.5 * dx, L + nghost * dx - 0.5 * dx, N)
for j in range(N):
    U[j] = InitialCondition(x[j])

# Visualization setup
fig, ax = plt.subplots(3, 1, sharex=True, sharey=False, dpi=140)
fig.subplots_adjust(hspace=0.1, wspace=0.0)

line_d, = ax[0].plot([], [], 'r-o', ls='-', markeredgecolor='k', markersize=3)
line_u, = ax[1].plot([], [], 'b-o', ls='-', markeredgecolor='k', markersize=3)
line_p, = ax[2].plot([], [], 'g-o', ls='-', markeredgecolor='k', markersize=3)

ax[2].set_xlabel('x')
ax[0].set_ylabel('Density')
ax[1].set_ylabel('Velocity')
ax[2].set_ylabel('Pressure')

ax[0].set_xlim(0.0, (L-0.25)*(1.0/0.75))
ax[0].set_ylim(+0.0, 1.5e3)
ax[1].set_ylim(-0.2, 1.0)
ax[2].set_ylim(+0.0, 700)

nframe = 99999999  # arbitrarily large
nstep_per_image = 1  # number of steps per image
anim = animation.FuncAnimation(fig, func=update, init_func=init, frames=nframe, interval=200, repeat=False)

plt.show()
