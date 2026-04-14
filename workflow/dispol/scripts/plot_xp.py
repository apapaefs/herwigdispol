import numpy as np
import plotext as plt


rho = np.linspace(0, 1,1000000)

n=0.2
xp2 = 1 - rho**(1/(1-n))
n=0.4
xp4 = 1 - rho**(1/(1-n))
n=0.6
xp6 = 1 - rho**(1/(1-n))
n=0.8
xp8 = 1 - rho**(1/(1-n))

plt.plot(rho, xp2)
plt.plot(rho, xp4)
plt.plot(rho, xp6)
plt.plot(rho, xp8)
plt.title("xp(rho)")
plt.show()
