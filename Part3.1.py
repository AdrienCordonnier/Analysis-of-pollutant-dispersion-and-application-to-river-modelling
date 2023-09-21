import numpy as np
import matplotlib.pyplot as plt

# Parameters
nu = V = 1
a = b = 50
T_min = 0
T_max = 5
h = 0.5
tau = 0.1
sigma = 1
N = int(a / h - h)
xe = ye = a / 2
T = T_max - T_min

def f(x, y):
    F = np.zeros(N ** 2)
    for i in range(N):
        for j in range(N):
            F[i * N + j] = (sigma * np.sqrt(2 * np.pi)) * np.exp(-((x[i] - xe) ** 2 + (y[j] - ye) ** 2)/ (2 * sigma ** 2))
    return F

def dispertion2D(N):
    x = np.linspace(h, a - h, N)
    y = np.linspace(h, a - h, N)
    u0 = f(x, y)

    sd = 1 / (h ** 2) * np.ones(N**2 - 1)
    d = -4 / h ** 2 * np.ones(N**2)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = -4 / h ** 2 * np.ones(N**2 - 1)
    d = 4 / h ** 2 * np.ones(N**2)
    M = np.diag(sd, -1) + np.diag(d, 0) + np.diag(-sd, 1)

    M1 = np.eye(N**2) - (nu * tau / 2) * A + nu * tau * M
    M2 = np.eye(N**2) + (nu * tau / 2) * A - nu * tau * M

    nmax = int(T / tau)
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = np.linalg.solve(M1, M2 @ u)
        l_u.append(u)

    return x, y, l_u


# Call the dispersion2D function
x, y, l_u = dispertion2D(N)
X, Y = np.meshgrid(x, y)

for i in l_u:
    plt.scatter(X, Y, c=i, cmap='jet')
    plt.colorbar(label='Intensité')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Graphe 2D')
    plt.show()
