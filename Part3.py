import numpy as np
import matplotlib.pyplot as plt


def build_matrix_M(N, h):
    M = np.zeros((N * N, N * N))
    for i in range(N * N):
        M[i, i] = 4 / h ** 2
        if i % N != 0:
            M[i, i - 1] = -1 / h ** 2
        if (i + 1) % N != 0:
            M[i, i + 1] = -1 / h ** 2
        if i >= N:
            M[i, i - N] = -1 / h ** 2
        if i < N * (N - 1):
            M[i, i + N] = -1 / h ** 2
    return M


def solve_laplace_eq(N, n, k, a, b):
    h = a / N

    xi = np.linspace(h, a, N)
    yj = np.linspace(h, b, N)

    M = build_matrix_M(N, h)

    F = np.zeros(N ** 2)

    for i in range(N):
        for j in range(N):
            x = xi[i]
            y = yj[j]
            F[i * N + j] = np.sin(n * np.pi * x / a) * np.sin(k * np.pi * y / b) * (
                        (n * np.pi / a) ** 2 + ((k * np.pi / b) ** 2))

    # Résolution du système linéaire MU = F
    U = np.linalg.solve(M, F)

    # Calcul de la solution exacte
    uex = np.zeros(N ** 2)
    for i in range(N):
        for j in range(N):
            x = xi[i]
            y = yj[j]
            uex[i * N + j] = np.sin(n * np.pi * x / a) * np.sin(k * np.pi * y / b)

    return U, uex


# Paramètres du domaine
a = b = 2

# Paramètres de la solution
n = 1
k = 2

# Paramètre de discrétisation
N = 160

# Résolution numérique
U, uex = solve_laplace_eq(N, n, k, a, b)

# Affichage des résultats
x = np.linspace(0, a, N)
y = np.linspace(0, b, N)
X, Y = np.meshgrid(x, y)

plt.figure(1)
plt.imshow(uex.reshape((N, N)), cmap='jet', extent=[0, a, 0, b])
plt.colorbar(label='Solution exacte')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution exacte de l\'équation de Laplace')

plt.figure(2)
plt.imshow(U.reshape((N, N)), cmap='jet', extent=[0, a, 0, b])
plt.colorbar(label='Solution numérique')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution numérique de l\'équation de Laplace')

### Etude de la convergence
# Liste pour stocker les erreurs
errors = []
# Valeurs de discrétisation à tester
N_values = [10, 20, 40, 80, 160]
# Calcul de l'erreur pour chaque valeur de N
for N in N_values:
    print(N)
    # Résolution de l'équation de Laplace
    U, uex = solve_laplace_eq(N, n, k, a, b)

    # Calcul de l'erreur d'approximation
    norm_error = np.linalg.norm(U.reshape((N, N)) - uex.reshape((N, N))) / N
    errors.append(norm_error)

# Affichage des résultats
plt.figure(3)
plt.plot(N_values, errors)
plt.xlabel('Nombre de points de discrétisation (N)')
plt.ylabel("Erreur d'approximation")
plt.title("Convergence de l'erreur d'approximation")
plt.grid()
plt.show()


