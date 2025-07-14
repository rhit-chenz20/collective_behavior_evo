import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create meshgrid of x and y values
x = np.linspace(-1, 1, 12)
y = np.linspace(-1, 1, 12)
X, Y = np.meshgrid(x, y)

# Define Z based on some function of X and Y
Z = (1 / (1 - X * Y))**2  # for example, like in Moore et al.

# Plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
ax.view_init(elev=30, azim=50)
# ax.view_init(elev=-2, azim=-90)

# Add labels
ax.set_xlabel('ψ12')
ax.set_ylabel('ψ21')
ax.set_zlabel('Relative change in phenotypic evolution')
ax.set_title('Relative Rate of Evolution (Two-Trait Reciprocal)')

plt.colorbar(surf, shrink=0.5, aspect=5)
plt.show()