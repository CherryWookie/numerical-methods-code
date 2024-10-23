import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

# Given parameters
a = 5.0        # Width of potential well
V0 = 50.0      # Potential inside the well
x0 = -20.0     # Initial position of the wave packet
sigma = 1.0    # Width of the initial wave packet
k0 = 10.0      # Wave number
L = 100.0      # Spatial domain [-L, L]
dx = 5e-2      # Spatial step size
dt = 5e-3      # Time step size
t_final = 5.0  # Final time
hbar = 1.0     # Reduced Planck constant (scaled to 1)
m = 1.0        # Mass of the particle (scaled to 1)

# Derived parameters
N = int(2 * L / dx) + 1  # Number of spatial points
x = np.linspace(-L, L, N)  # Spatial grid
steps = int(t_final / dt)  # Number of time steps

# Potential function V(x)
def potential(x, V0, a):
    V = np.zeros_like(x)
    V[np.abs(x) <= a / 2] = V0
    return V

V = potential(x, V0, a)

# Initial wave packet Psi(0, x)
def psi_initial(x, x0, k0, sigma):
    normalization = 1 / (np.sqrt(np.sqrt(2 * np.pi) * sigma))
    psi_0 = normalization * np.exp(- (x - x0)**2 / (4 * sigma**2)) * np.exp(1j * k0 * (x - x0))
    return psi_0

psi = psi_initial(x, x0, k0, sigma)

# Crank-Nicolson coefficients
r = dt / (2 * dx**2)  # r = dt / (2 * dx^2)

# Hamiltonian matrix H (tridiagonal matrix)
H = np.zeros((N, N), dtype=complex)

for i in range(1, N - 1):
    H[i, i] = -2 * r + 1j * dt * V[i] / (2 * hbar)
    H[i, i - 1] = r
    H[i, i + 1] = r

# Impose boundary conditions: Psi(-L) = Psi(L) = 0
H[0, 0] = H[N - 1, N - 1] = 1
H[0, 1] = H[N - 1, N - 2] = 0

# Crank-Nicolson matrices
A = np.eye(N, dtype=complex) - 1j * H  # A = I - i * H
B = np.eye(N, dtype=complex) + 1j * H  # B = I + i * H

# Time evolution using Crank-Nicolson
def crank_nicolson_step(psi, A, B):
    # Solve A Psi^{n+1} = B Psi^n
    rhs = B @ psi
    psi_new = la.solve(A, rhs)
    return psi_new

# Time evolution loop
psi_t = [psi.copy()]  # To store Psi at each time step for plotting

for step in range(steps):
    psi = crank_nicolson_step(psi, A, B)
    psi_t.append(psi.copy())

# Convert psi_t to a numpy array for easier handling
psi_t = np.array(psi_t)

# Visualization: Initial and final wave function
def plot_wave_function(psi, x, step):
    plt.plot(x, np.abs(psi)**2, label=f"t = {step * dt:.2f}")
    plt.xlabel("x")
    plt.ylabel("|Psi(x)|^2")
    plt.legend()
    plt.grid()

plt.figure(figsize=(8, 5))
plot_wave_function(psi_t[0], x, 0)
plot_wave_function(psi_t[-1], x, steps)
plt.title("Wave Function at Initial and Final Time Steps")
plt.show()

# Optional: create an animation of the wave packet evolution
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
line, = ax.plot(x, np.abs(psi_t[0])**2)

def update(frame):
    line.set_ydata(np.abs(psi_t[frame])**2)
    ax.set_title(f"Wave Function at t = {frame * dt:.2f}")
    return line,

ani = FuncAnimation(fig, update, frames=range(0, steps, 10), blit=True)
plt.show()
