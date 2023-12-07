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

    return [dCdt, dPdt, dRdt, dTdt]

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

lambda_c_values = np.array([ 10**(-6), 10**(-7), 10**(-8), 10**(-9)])

# Time vector
t0 = 0
T = 500
pas = 0.001
t = np.arange(t0,T,pas)

fig, axs = plt.subplots(2, 2, figsize=(12, 8))

for i, ax in enumerate(axs.flatten()):

    lambda_c = lambda_c_values[i]
    # Solving the EDO system
    untreated = spi.odeint(model, x0, t, args=(mu_c, mu_p, gamma_p, gamma_c, k_t))
    TGF_beta = spi.odeint(model, x0, t, args=(mu_c, 0.9*mu_p, 0.9*gamma_p, 0.9*gamma_c, k_t))
    immune_act = spi.odeint(model, x0, t, args=(mu_c, mu_p, gamma_p, gamma_c, 2*k_t))
    both_s = spi.odeint(model, x0, t, args=(mu_c, 0.9*mu_p, 0.9*gamma_p, 0.9*gamma_c, 2*k_t))

    # Graph
    ax.plot(t, untreated[:, 0], label='Untreated patient')
    ax.plot(t, TGF_beta[:, 0], label='TGF_beta silencing')
    ax.plot(t, immune_act[:, 0], label='Immune activation')
    ax.plot(t, both_s[:, 0], label='Both')
    
    ax.set_title('λc = ' + str(lambda_c))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Cancer cells (C)')
    ax.legend()
plt.tight_layout()
plt.show()

#Plotting graph for lambda_c = 1e-7 with all 4 solutions of the model
lambda_c = 1e-7
plt.figure()
plt.clf()

plt.subplot(221)
plt.plot(t,untreated[:,0],'.-')
plt.ylabel('PCC')
plt.xlabel('t')
plt.subplot(222)
plt.plot(t,untreated[:,1])
plt.ylabel('PsC')
plt.xlabel('t')
plt.subplot(223)
plt.plot(t,untreated[:,2])
plt.ylabel('R')
plt.xlabel('t')
plt.subplot(224)
plt.plot(t,untreated[:,3])
plt.ylabel('T') 
plt.xlabel('t')
    
plt.legend('Untreatd')    
plt.suptitle('lambda_c='+str(lambda_c))
plt.show()