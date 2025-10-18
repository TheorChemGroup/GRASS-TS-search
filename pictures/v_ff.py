# Import required modules 
import numpy as np 
import matplotlib.pyplot as plt
e=2.718281828 
#1f=lambda x,y :-x**2 + y**2
#2 f=lambda x,y :-(x+0.5*y)**2 + y**2
f=lambda x,y :np.cos(3*0.31415*x) + np.cos(3*0.31415*y)
d_x=lambda x,y,f:(f(x,y)-f(x+0.001,y))/0.001
d_y=lambda x,y,f:(f(x,y)-f(x,y+0.001))/0.001
# Meshgrid 
x, y = np.meshgrid(np.linspace(-5, 5, 10),  
                   np.linspace(-5, 5, 10)) 
xc, yc = np.meshgrid(np.linspace(-6, 6, 100),  
                   np.linspace(-6, 6, 100)) 

# Directional vectors 
u = d_x(x,y,f)
v = -d_y(x,y,f)

plt.axes().set_aspect(1)  
plt.grid(zorder=0)

# Plotting Vector Field with QUIVER 
plt.quiver(x, y, u, v, color='g',zorder=6) 
#plt.contour(xc, yc, f(xc,yc),zorder=5) 
plt.title('') 
  
#plt.scatter([0,0], [3.333,-3.333], color="r",zorder=7)
plt.scatter([3.333,-3.333], [0,0], color="orange",zorder=7)

# Setting x, y boundary limits 
plt.xlim(-6.5, 6.5) 
plt.ylim(-6.5, 6.5) 
  
# Show plot with grid 
 
plt.savefig("pictures/sadd3_rev_y_ff",dpi=300)