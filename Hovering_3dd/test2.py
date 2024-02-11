import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
from IPython.display import display, clear_output

t = np.linspace(0, 2*np.pi)
x = np.sin(t)

fig, ax = plt.subplots()
l, = ax.plot([0, 2*np.pi], [-1, 1])

animate = lambda i: l.set_data(t[:i], x[:i])

for i in range(len(x)):
    animate(i)
    clear_output(wait=True)
    display(fig)

plt.show()
