import matplotlib.pyplot as plt  # Tracer des graphes
import numpy as np               # Calculer/générer des matrices
import imageio                   # Générer des GIFs
import os                        # Spécifier les différents chemins pour l'enregistrement

sigma = 1
xe = 25

#Modifiez ce chemin pour enregistrer les gifs
chemin_GIF = "C:/Users/musta/Documents/IPSA-Cours/Aero3/Semestre2/Ma323 - Méthode des différences finies/Projet/GIF"

def f(x):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x-xe)**2/(2*sigma**2))


def Solexplicitecentre():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = V / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    Me = np.eye(N) - tau * M

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = Me @ u
        l_u.append(u)
    return x, l_u


def SolCN():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M
    M2 = np.eye(N) - V * tau / 2 * M
    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = np.linalg.solve(M1, M2 @ u)
        l_u.append(u)

    return x, l_u


def SolCN2():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = 1 / h**2 * np.ones(N - 1)
    d = -2/h**2 * np.ones(N)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M - nu * tau / 2 * A
    M2 = np.eye(N) - V * tau / 2 * M + nu * tau / 2 * A

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = np.linalg.solve(M1, M2 @ u)
        l_u.append(u)

    return x, l_u


def SolCN3():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = 1 / h**2 * np.ones(N - 1)
    d = -2/h**2 * np.ones(N)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M - nu * tau / 2 * A
    M2 = np.eye(N) - V * tau / 2 * M + nu * tau / 2 * A

    M1[0][N-1] = M1[1][0]
    M1[N-1][0] = M1[0][1]
    M2[0][N - 1] = M2[1][0]
    M2[N - 1][0] = M2[0][1]

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = np.linalg.solve(M1, M2 @ u)
        l_u.append(u)

    return x, l_u


def SolCN4():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = 1 / h**2 * np.ones(N - 1)
    d = -2/h**2 * np.ones(N)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M - nu * tau / 2 * A
    M2 = np.eye(N) - V * tau / 2 * M + nu * tau / 2 * A

    M1[0][N-1] = M1[1][0]
    M1[N-1][0] = M1[0][1]
    M2[0][N - 1] = M2[1][0]
    M2[N - 1][0] = M2[0][1]

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        u = np.linalg.solve(M1, M2 @ u + tau * u0) # Ajout du terme source
        l_u.append(u)

    return x, l_u


def SolCN5():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = 1 / h**2 * np.ones(N - 1)
    d = -2/h**2 * np.ones(N)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M - nu * tau / 2 * A
    M2 = np.eye(N) - V * tau / 2 * M + nu * tau / 2 * A

    M1[0][N-1] = M1[1][0]
    M1[N-1][0] = M1[0][1]
    M2[0][N - 1] = M2[1][0]
    M2[N - 1][0] = M2[0][1]

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    for n in range(nmax):
        if n % 2 == 0:
            u = np.linalg.solve(M1, M2 @ u + u0 * tau)  # Ajout du terme source
        else:
            u = np.linalg.solve(M1, M2 @ u)
        l_u.append(u)

    return x, l_u


def SolCN6():
    x = np.linspace(h, a-h, N)
    u0 = f(x)

    sd = 1 / (2*h) * np.ones(N - 1)
    d = np.zeros(N)
    M = np.diag(-sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    sd = 1 / h**2 * np.ones(N - 1)
    d = -2/h**2 * np.ones(N)
    A = np.diag(sd, -1) + np.diag(d, 0) + np.diag(sd, 1)

    M1 = np.eye(N) + V * tau / 2 * M - nu * tau / 2 * A
    M2 = np.eye(N) - V * tau / 2 * M + nu * tau / 2 * A

    M1[0][N-1] = M1[1][0]
    M1[N-1][0] = M1[0][1]
    M2[0][N - 1] = M2[1][0]
    M2[N - 1][0] = M2[0][1]

    nmax = int(T / tau)  # t = n*tau d'où n = t/tau
    u = u0
    l_u = [u]
    fonctionnement, pause = 1, 0
    for n in range(nmax):
        u1 = np.linalg.solve(M1, M2 @ u + u0 * tau)
        if np.max(u1) <= 4/10:
            u = u1
            fonctionnement += 1
        else:
            u = np.linalg.solve(M1, M2 @ u)
            pause += 1
        l_u.append(u)

    return x, l_u, fonctionnement * tau, pause * tau

def gif(x, l_u, methode, cas):
    if cas == "a":
        d=0.1               # Vitesse du GIF
        step = 50           # Le Pas pour les plots
    elif cas == "b":
        d=0.1
        step = 200
    else:
        d=0.1
        step = 300

    # Créer un nouveau dossier pour enregistrer les images
    if not os.path.exists("Graphe_{}_cas_{}".format(methode, cas)):
        os.makedirs("Graphe_{}_cas_{}".format(methode, cas))
        os.chdir("Graphe_{}_cas_{}".format(methode, cas))
    else:
        os.makedirs("Graphe_{}_cas_{} - Copie".format(methode, cas))
        os.chdir("Graphe_{}_cas_{} - Copie".format(methode, cas))

    # Enregistrer la figure initiale
    fig, ax = plt.subplots()
    plt.figure(figsize=(12, 8))
    ax.plot(x, l_u[0])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    tx1, tx2 = "$\dfrac{V\u03C4}{h}$", "$\dfrac{\u03BD\u03C4}{h}$"
    ax.set_title(f"Résolution à l'aide de : {methode}\n$h = {h} ; \u03C4 = {tau}$ ; {tx1} = {c_V} ; {tx2} = {c_nu}")
    ax.grid(True, color='gray', linestyle='--', linewidth=0.5)
    fig.savefig('fig_0.png')

    # Enregistrer les autres figures
    figs = []
    print("Plotting...")
    for i in range(2, len(l_u), step):
        fig, ax = plt.subplots()
        plt.figure(figsize=(12, 8))
        ax.plot(x, l_u[0], label='t = 0')
        ax.plot(x, l_u[i], label=f't = {i}')
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.legend()
        tx1, tx2 = "$\dfrac{V\u03C4}{h}$", "$\dfrac{\u03BD\u03C4}{h}$"
        ax.set_title(f"Résolution à l'aide de : {methode}\n$h = {h} ; \u03C4 = {tau}$ ; {tx1} = {c_V} ; {tx2} = {c_nu}")
        ax.grid(True, color='gray', linestyle='--', linewidth=0.5)
        fig.savefig(f'fig_{i - 1}.png')
        figs.append(imageio.imread(f'fig_{i - 1}.png'))

    # Créer le GIF à partir des images enregistrées
    imageio.mimsave('animation.gif', figs, duration=d, loop=500)
    os.rename('animation.gif', chemin_GIF + "/Graphe_{}_cas_{}.gif".format(methode, cas))
    os.chdir("..")


all_methode = {"Explicite Centré":Solexplicitecentre,
               "Crank-Nicolson": SolCN,
               "Crank-Nicolson - Convection diffusion": SolCN2,
               "Crank-Nicolson (Neumann)": SolCN3,
               "Crank-Nicolson avec terme source" : SolCN4,
               "Crank-Nicolson avec terme source que le jour" : SolCN5,
               "Crank-Nicolson avec régulation du polluant" : SolCN6}


Cas = {"a":{"a":50, "T":5, "V":1, "nu":1, "h":0.1, "tau":0.0025},
       "b":{"a":50, "T":25, "V":1, "nu":1, "h":0.1, "tau":0.0025},
       "c":{"a":100, "T":75, "V":1, "nu":1, "h":0.1, "tau":0.0025}}


for cas in Cas:
    if cas == "a" or cas == "b":
        for methode in all_methode:
            a = Cas[cas]["a"]
            T = Cas[cas]["T"]
            V = Cas[cas]["V"]
            nu = Cas[cas]["nu"]
            h = Cas[cas]["h"]
            tau = Cas[cas]["tau"]

            N = int((a-h)/h)
            c_V = round(V * tau / h, 3)
            c_nu = round(nu * tau / h, 3)

            if methode == "Crank-Nicolson avec régulation du polluant" :
                x, l_u, fonctionnement, pause = all_methode[methode]()
                gif(x, l_u, methode, cas)
                print(f"temps de fonctionnement : {fonctionnement}\ntemps de pause : {pause}")
            else :
                x, l_u = all_methode[methode]()
                gif(x, l_u, methode, cas)

    if cas == "c":
        a = Cas[cas]["a"]
        T = Cas[cas]["T"]
        V = Cas[cas]["V"]
        nu = Cas[cas]["nu"]
        h = Cas[cas]["h"]
        tau = Cas[cas]["tau"]

        N = int((a - h) / h)
        c_V = round(V * tau / h, 3)
        c_nu = round(nu * tau / h, 3)
        
        x, l_u, fonctionnement, pause = SolCN6()
        gif(x, l_u, "Crank-Nicolson avec régulation du polluant", cas)
        print(f"temps de fonctionnement : {fonctionnement}\ntemps de pause : {pause}")









