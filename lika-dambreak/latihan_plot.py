import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

x = np.linspace(0,24*np.pi,512)
y = np.sin(x)

def fft(x):
    fft = np.abs(np.fft.rfft(x))
    return fft**2/(fft**2).max()

fig, (ax1,ax2) = plt.subplots(nrows=2)
line1, = ax1.plot(x,y)
line2, = ax2.plot(fft(y))
ax2.set_xlim(0,50)
ax2.set_ylim(0,1)

def update(i):
    y = np.sin((i+1)/30.*x)
    line1.set_data(x,y)
    y2 = fft(y)
    line2.set_data(range(len(y2)), y2)

ani = matplotlib.animation.FuncAnimation(fig, update, frames=60, repeat=True)
plt.show()
