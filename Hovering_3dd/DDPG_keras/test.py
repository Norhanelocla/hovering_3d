import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3D
import numpy as np
import csv
import time

with open('drone_x.csv', 'r') as f:
    reader = csv.reader(f)
    data = list(reader)
drone_x = np.array(data, dtype=float)

with open('drone_y.csv', 'r') as f:
    reader = csv.reader(f)
    data2 = list(reader)
drone_y = np.array(data2, dtype=float)

with open('drone_z.csv', 'r') as f:
    reader = csv.reader(f)
    data3 = list(reader)
drone_z = np.array(data3, dtype=float)

print(drone_x)
window_x = [0 ,0 ,0 ,0]
window_y = [1 ,2 ,3 ,4]
window_z = [0 ,0 ,0 ,0]

class YourClassNameHere:
    def __init__(self,drone_x,drone_y,drone_z,window_x,window_y,window_z):
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.counter = 0
        self.drone_x = drone_x
        self.drone_y = drone_y
        self.drone_z = drone_z
        self.window_x = window_x
        self.window_y = window_y
        self.window_z = window_z

    def animate(self, num_frames):
        for frame in range(num_frames):
            self.counter += 1
            #print(self.window_x[frame])
            #self.ax.cla()

            # Plot the window's 3D position
            drone_position = self.ax.scatter(self.drone_x[frame], self.drone_y[frame], self.drone_z[frame], c='r', marker='o', label='Window Position')

            # Set labels and limits for the 3D plot
            self.ax.set_xlabel('X (m)')
            self.ax.set_ylabel('Y (m)')
            self.ax.set_zlabel('Z (m)')
            self.ax.set_xlim([-100, 100])
            self.ax.set_ylim([-100, 100])
            self.ax.set_zlim([-100, 100])
            self.ax.legend()

            # Redraw the 3D plot
            plt.draw()
            plt.pause(0.1)

    def render(self, mode='human', close=False):
        pass  # Leave this empty for animation


animation_example = YourClassNameHere(drone_x,drone_y,drone_z,window_x,window_y,window_z)
# animation_example.animate(num_frames=len(time_vals))
animation_example.animate(num_frames=1)
plt.show()
