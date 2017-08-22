import numpy as np
import matplotlib.pyplot as plt

xbins, freq, x, prob = np.genfromtxt("hist.dat").T
Ns, means, means2 = np.genfromtxt("N.dat").T

m, b = np.polyfit(Ns, means2, 1)

dx = x[1]-x[0]

fig, ax1 = plt.subplots()
ax1.bar(xbins, freq, 2.5, color='b')

ax1.set_xlabel('$x$')
ax1.tick_params('y', colors='b')
ax1.set_ylabel('Freq', color='b')

ax2 = ax1.twinx()
ax2.plot(x, prob, c="r")
ax2.tick_params('y', colors='r')
ax2.set_ylabel('Probability', color='r')
ax2.set_ylim(0, max(prob))

fig.tight_layout()
fig.savefig("hist.png")

fig, ax1 = plt.subplots()

x = np.logspace(0, 1.1*np.log10(Ns[-1]), 10)
y = m*x + b

ax1.loglog(Ns, means2, color = 'b')
ax1.loglog(x, y, "--", color = 'k', lw=1)
ax1.text(5e1, 10e3, "$<x^2> = %.3fN + %.3f$"%(m, b))

ax1.tick_params('y', colors='b')
ax1.set_ylabel('$<x^2>$', color='b')

ax2 = ax1.twinx()
ax2.plot(Ns, means, color = 'r')
ax2.tick_params('y', colors='r')
ax2.set_ylabel('<x>', color='r')
ax1.set_xlabel('$N$ steps')

fig.tight_layout()
fig.savefig("means.png")
