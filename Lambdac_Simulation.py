import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

def model(x, t, mu_c, mu_p, gamma_p, gamma_c, k_t):

    C, P, R, T = x

    # Equations arguments
    arg11 = (k_c + (mu_c * P)) * C**(3/4) * (1 - (C/C0)**(1/4))
    arg12 = (lambda_c * C * T) / (K_c + (1 - R))
    arg21 = (k_p + ((mu_p * C) / (K_p + C))) * P * (1 - P/P0)
    arg22 = lambda_p * P
    arg3 = k_r - ((lambda_r + (gamma_p * P) + (gamma_c * C)) * R)
    arg41 = (k_t * R) / (K_t + (1 - R))
    arg42 = lambda_t * T

    # EDOs
    dCdt = arg11 - arg12
    dPdt = arg21 - arg22
    dRdt = arg3
    dTdt = arg41 - arg42

    return dCdt, dPdt, dRdt, dTdt

# Initial conditions
x0 = [10,100,0.5,300]

# Parameters of the system
C0=10**6
P0=10**5
k_c=0.075
mu_c=20*k_c/P0
K_c=0.1
k_p=0.2
mu_p=20*k_p
K_p=C0/100
lambda_p=0.15
k_r=0.2
lambda_r=0.22
Ps=P0*(1-lambda_p/k_p)
gamma_p=0.02*lambda_r/Ps
gamma_c=gamma_p
k_t=3300
K_t=K_c
lambda_t=0.3

# Time vector
t0 = 0
T = 300
pas = 1
t = np.arange(t0, T, pas)

# Values of lambda
lambda_values = np.logspace(-9, -6.5, 50)

# Store results
dCdt_values_un = []
dCdt_values_TGF = []
dCdt_values_immune = []
dCdt_values_both = []
dCdt_values_article = []


# Run the model for different values of lambda
for lambda_c in lambda_values:
    
    untreated = spi.odeint(model, x0, t, args=(mu_c, mu_p, gamma_p, gamma_c, k_t))
    article = spi.odeint(model, x0, t, args=(mu_c, mu_p, 0*gamma_p, 0*gamma_c, k_t))
    
    dCdt_values_un.append(untreated[-1, 0])  # Store the value of dC/dt at the final time point
    dCdt_values_article.append(article[-1, 0])

# Plot the graph
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(lambda_values, dCdt_values_un, label = "Untreated")
ax.plot(lambda_values, dCdt_values_article, label = "gama_p = gama_c = 0")

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('λc')
ax.set_ylabel('Cancer cells (C)')
ax.legend()
plt.show()
