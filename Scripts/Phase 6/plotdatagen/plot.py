import matplotlib.pyplot as plt
import numpy as np
plt.style.use("seaborn")
data = np.loadtxt("data(const).dat")

x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]
y4 = data[:, 4]
y5 = data[:, 5]
y6 = data[:, 6]
y7 = data[:, 7]
y8 = data[:, 8]
y9 = data[:, 9]
y10 = data[:, 10]
y11 = data[:, 11]

plt.xlabel("scaling constant of the GTOs")
plt.ylabel("vector entry magnitude")
plt.plot(x, y1, label="1")
plt.plot(x, y2, label="2")
plt.plot(x, y3, label="3")
plt.plot(x, y4, label="4")
plt.plot(x, y5, label="5")
plt.plot(x, y6, label="6")
plt.plot(x, y7, label="16")
plt.plot(x, y8, label="19")
plt.plot(x, y9, label="20")
plt.plot(x, y10, label="25")
plt.plot(x, y11, label="30")
plt.vlines(0.19, -0.1, 5, linestyles='dashed', color='red', label=r"const=0.19")
plt.legend(loc="upper right")
plt.title("FP vector entries vs. scaling constant of the GTOs")
plt.xlim(0, 2.5)
plt.ylim(-0.1, 5)
plt.savefig("plot3.png")