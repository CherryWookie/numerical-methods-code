import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg



c = 100

x = np.linspace(0, 1, 101)
a = 1/100 # Spacial Grid

time = np.linspace(0, 1, 1001)
h = 0.01/1000 # Spacing (delta t I think)

y = np.zeros(101) # Horizontal, fixed at end points
v = x*(1-x)*np.exp(-(x-0.5)**2/100) # Pertubation? Striking the wave... might be unecessary 

# Set up banded matrix
k = 0.25*h**2*c**2/a**2 # Constant with hc/4a or whatever
A1 = 1 + 2*k
A2 = -k
A = np.zeros((3, 99)) # exclude boundary conditions. 99 = 101 - 2
A[0] = A2
A[1] = A1
A[2] = A2
 

shape = []
for t in time:
    shape.append(y.copy())
    y_old = y.copy() # Copy of old wavefunction
    b = k*y[:-2] + (1-2*k)*y[1:-1] + k*y[2:] + h*v[1:-1]
    y[1:-1] = linalg.solve_banded((1,1), A, b) # (1,1) = one line above and one line below diagonal. A = array, b value
    v[1:-1] += 2*k/h*(y[:-2] + y[2:] -2*y[1:-1] + y_old[:-2] + y_old[2:] - 2*y_old[1:-1])
    
shape = np.array(shape)
plt.plot((time, shape[:, 51])) # middle element
plt.show()
