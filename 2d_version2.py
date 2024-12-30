import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('out.txt', skiprows=1) 
X = data[1:, 0]  
Y = data[1:, 1]  
Z = data[1:, 3] 

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X, Y, Z, color='blue', s=50, label='Data points')

ax.set_title('3D Plot ', fontsize=16)
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.legend()


plt.show()