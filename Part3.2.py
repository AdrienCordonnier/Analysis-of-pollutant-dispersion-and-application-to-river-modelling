import numpy as np
import matplotlib.pyplot as plt

# Parameters
V = 1
a = b = 50
T_min = 0
T_max = 5
h = 0.5
tau = 0.1
sigma = 1

# Grid dimensions
N_x = int((a - 0) / h) + 1
N_y = int((b - 0) / h) + 1
N_t = int((T_max - T_min) / tau) + 1

# Matrices to store data
sd = V / (2 * h) * np.ones(N - 1)
d = np.zeros(N)
M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    Me = np.eye(N) - tau * M
Mx = np.zeros((N_x * N_y, N_t))
My = np.zeros((N_x * N_y, N_t))

# Matrix A
A = np.zeros((N_x * N_y, N_x * N_y))
for i in range(N_x * N_y):
    A[i, i] = -4 / h**2
    if i - 1 >= 0 and (i % N_y) != 0:
        A[i, i - 1] = -1 / h**2
    if i + 1 < N_x * N_y and (i % N_y) != (N_y - 1):
        A[i, i + 1] = -1 / h**2
    if i - N_y >= 0:
        A[i, i - N_y] = -1 / h**2
    if i + N_y < N_x * N_y:
        A[i, i + N_y] = -1 / h**2

# Initial conditions
x_e = y_e = a / 2
x = np.linspace(0, a, N_x)
y = np.linspace(0, b, N_y)
X, Y = np.meshgrid(x, y)
Mx[:, 0] = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((X - x_e) ** 2 + (Y - y_e) ** 2) / (2 * sigma ** 2)).flatten()
My[:, 0] = 0

# Numerical simulation
for k in range(1, N_t):
    Mx[:, k] = Mx[:, k - 1] + tau * (V / h ** 2) * (np.dot(A, Mx[:, k - 1]) - 4 * Mx[:, k - 1])
    My[:, k] = My[:, k - 1] + tau * (V / h ** 2) * (np.dot(A, My[:, k - 1]) - 4 * My[:, k - 1])

print(A)

# Approximation of the continuous solution
M = np.sqrt(Mx**2 + My**2)/10**55

# Visualization of the results
fig, ax = plt.subplots()
contour = ax.contourf(X, Y, M[:, -1].reshape(N_x, N_y), cmap='jet')
cbar = fig.colorbar(contour, ax=ax, label='Concentration')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Dispersion du polluant dans eau')
ax.set_xlim([a/4, 3*a/4])
ax.set_ylim([b/4, 3*b/4])
plt.show()


