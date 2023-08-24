import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters
beta = 0.3       # Constant transmission rate
gamma = 0.1      # Recovery rate
lambda_val = 0.05  # Immunity loss rate
alpha = 0.08     # Immunity waning rate
sigma = 0.005    # Vaccination rate (fraction of susceptibles vaccinated per day)
delta = 0.01     # Demographic balanced birth and death rate
initial_population = 1000  # Total population size

# System of ODEs
def sirs_model(y, t, beta, gamma, lambda_val, alpha, sigma, delta):
    S, I, R, V = y

    dSdt = delta - (beta * S * I) / (S + I + R + V) + lambda_val * R - sigma * S
    dIdt = (beta * S * I) / (S + I + R + V) - gamma * I - alpha * I - delta * I
    dRdt = gamma * I + alpha * I - lambda_val * R - delta * R
    dVdt = sigma * S + alpha * R - delta * V

    return [dSdt, dIdt, dRdt, dVdt]

# Initial conditions
S0 = 0.99 * initial_population
I0 = 0.01 * initial_population
R0 = 0.0
V0 = 0.0

# Time points
t = np.linspace(0, 100, 1000)

# Solve the ODEs
solution = odeint(sirs_model, [S0, I0, R0, V0], t, args=(beta, gamma, lambda_val, alpha, sigma, delta))
S, I, R, V = solution.T

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, S, label='Susceptible')
plt.plot(t, I, label='Infectious')
plt.plot(t, R, label='Recovered')
plt.plot(t, V, label='Vaccinated')
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('SIRS Model Simulation as ODEs')
plt.legend()
plt.grid()
plt.show()
