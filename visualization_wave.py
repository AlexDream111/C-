import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("wave_solution.dat")
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]
nx = int(np.sqrt(len(x)))

X, Y = x.reshape(nx, nx), y.reshape(nx, nx)
U = u.reshape(nx, nx)

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, U, levels=50, cmap='viridis')
plt.colorbar(label="Wave Amplitude")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Wave Equation Solution")
plt.savefig("wave_solution.png")
print("图像已保存为 wave_solution.png")



