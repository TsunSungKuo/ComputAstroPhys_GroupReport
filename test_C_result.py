import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 128 # real calculation zone
SS = 8
L = 1.0
nstep_per_image = 1     # plotting frequency
dx = L/N

x = np.linspace(0, L, N)
for j in range( N ):
   x[j] = (j+0.5)*dx    # cell-centered coordinates

# update animation

f = open('output.txt','r')
def init():
    f = open('output.txt','r')
    line_d.set_xdata( x )
    line_u.set_xdata( x )
    line_p.set_xdata( x )
   
   
    return line_d, line_u, line_p

def update( frame ):
    #f = open('output.txt','r')
    rho = np.zeros(N)
    vx  = np.zeros(N)
    P   = np.zeros(N)
    
    
    #for step in range(8):
    timeline = f.readline()    #read timestep
    for i in range(N):
        line = f.readline()
        #print(line)
        s = line.split(' ')
        rho[i] = s[0]
        vx[i]  = s[1]
        P[i]   = s[4]
        #print(s[4])
    line_d.set_ydata(rho)
    line_u.set_ydata(vx)
    line_p.set_ydata(P)
    ax[0].set_title(timeline)
    return line_d, line_u, line_p
            
t = 0.0

# create figure
fig, ax = plt.subplots( 3, 1, sharex=True, sharey=False, dpi=140 )
fig.subplots_adjust( hspace=0.1, wspace=0.0 )

line_d, = ax[0].plot( [], [], 'ro', markersize=2 )
line_u, = ax[1].plot( [], [], 'bo', markersize=2 )
line_p, = ax[2].plot( [], [], 'go', markersize=2 )
ax[2].set_xlabel( 'x' )
ax[0].set_ylabel( 'Density' )
ax[1].set_ylabel( 'Velocity' )
ax[2].set_ylabel( 'Pressure' )
ax[0].set_xlim( -0.1, L+0.1 )
ax[0].set_ylim( +0.0, 1.5e3 )
ax[1].set_ylim( -0.2, 1.0 )
ax[2].set_ylim( +0.0, 6.0e2 )


# plot theoretial result
file = 'strong_shock.txt'
data = np.loadtxt(file, dtype='float', delimiter=' ', skiprows=7, usecols=[0,2,4,10])
line_dt, = ax[0].plot( data[:,0], data[:,1], 'k-', ls='-', linewidth=1.0 )
line_ut, = ax[1].plot( data[:,0], data[:,2], 'k-', ls='-', linewidth=1.0 )
line_pt, = ax[2].plot( data[:,0], data[:,3], 'k-', ls='-', linewidth=1.0)

# create movie
nframe = 9999 # arbitrarily large

anim   = animation.FuncAnimation( fig, func=update, init_func=init,
                                  frames=nframe, interval=200, repeat=False )

plt.show()
f.close()

