import numpy as np
from scipy.linalg import solve_banded
from scipy import linalg


# Variables
a = 5
V_0 = 50 # Potential energy
x_0 = -20 # Initial position
sigma = 1
k_0 = 10 # Wavenumber
L = 100 # Length
dx = 5e-2 # Spatial mesh interval size
dt = 5e-3 # Timesteps
t_init = 0
t_final = 5
hbar = 1
m = 1
i = np.sqrt(-1)

# Generate Mesh
N = int(2 * L / dx) + 1  # Number of spatial points
x = np.linspace(-L, L, N)  # Spatial grid points
t_steps = int(t_final / dt)  # Number of time steps

# Initial Conditions for Psi(0,x)
norm_const = (np.sqrt(2*np.pi)**(1/2) * sigma) # Normalization Constant for ease of writing
psi_0 = 1 / ( norm_const * np.exp**((-(x - x_0)**2) / (4*sigma**2)) * np.exp**(i * k_0 *(x - x_0)))

# Potential

V = np.zeros(1,N+2)
idx1 = int((L - a / 2) / dx + 2)
idx2 = int((L + a / 2) / dx)

V[idx1:idx2 + 1] = V_0 # + 1 to include endpoint



# Coefficients
r = i * dt / 4 / dx**2
b = dt / dx**2

# First Matrix
a = (4 * i + (2 * dt * V[2:N + 1] / dx**2))
b = -dt / dx**2 * np.ones(1,N)
c = -dt / dx**2 * np.ones(1,N)

# LU solve

# Second Matrix
a = -a
b = -b
c = -c
