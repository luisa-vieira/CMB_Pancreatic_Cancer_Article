import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

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

# pParameters of the system
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
T = 500
pas = 0.001
t = np.arange(t0, T, pas)

# Number of simulations
num_simulations = 50

# Store results
untreated = []
TGF_beta = []
immune_act = []
both = []
llambda_c = 10**np.random.uniform(-9, -6.5, num_simulations)

# Run simulations for random lambda_c values
for i in range(num_simulations):

    lambda_c = llambda_c[i]

    untreated_si = spi.odeint(model, x0, t, args=(mu_c, mu_p, gamma_p, gamma_c, k_t))
    TGF_beta_si = spi.odeint(model, x0, t, args=(mu_c, 0.9*mu_p, 0.9*gamma_p, 0.9*gamma_c, k_t))
    immune_act_si = spi.odeint(model, x0, t, args=(mu_c, mu_p, gamma_p, gamma_c, 2*k_t))
    both_si = spi.odeint(model, x0, t, args=(mu_c, 0.9*mu_p, 0.9*gamma_p, 0.9*gamma_c, 2*k_t))

    untreated.append(untreated_si[:,0])
    TGF_beta.append(TGF_beta_si[:,0])
    immune_act.append(immune_act_si[:,0])
    both.append(both_si[:,0])


# Plot mean values
fig, ax = plt.subplots(figsize=(10, 6))

mean_untreated = np.mean(untreated, axis=0)
mean_TGF_beta = np.mean(TGF_beta, axis=0)
mean_immune_act = np.mean(immune_act, axis=0)
mean_both = np.mean(both, axis=0)

ax.plot(t, mean_untreated, label=f'Untreated')
ax.plot(t, mean_TGF_beta, label=f'TGF_beta silencing')
ax.plot(t, mean_immune_act, label=f'Immune activation')
ax.plot(t, mean_both, label=f'Both')
ax.set_title('Mean Values for Different λc')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Mean Cancer cells (C)')
ax.legend()
plt.show()

