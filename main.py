import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
beta = 0.3   # Transmission rate
gamma = 0.1  # Recovery rate
lambda_val = 0.05  # Immunity loss rate

# Initial conditions
S0 = 0.99  # Initial susceptible fraction
I0 = 0.01  # Initial infectious fraction
R0 = 0.0   # Initial recovered fraction

# Time points to simulate
t = np.linspace(0, 200, 1000)

# SIRS model equations
def sirs_model(y, t, beta, gamma, lambda_val):
    S, I, R = y
    dSdt = -beta * S * I + lambda_val * R
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I - lambda_val * R
    return [dSdt, dIdt, dRdt]

# Solve the differential equations
initial_conditions = [S0, I0, R0]
solution = odeint(sirs_model, initial_conditions, t, args=(beta, gamma, lambda_val))

S, I, R = solution.T

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, S, label='Susceptible')
plt.plot(t, I, label='Infectious')
plt.plot(t, R, label='Recovered')
plt.xlabel('Time')
plt.ylabel('Fraction of Population')
plt.title('SIRS Model Simulation')
plt.legend()
plt.grid()
plt.show()
