# Parameters
beta = 0.3  # Transmission rate
gamma = 0.1  # Recovery rate
lambda_val = 0.05  # Immunity loss rate

# Initial conditions
S = 0.99  # Initial susceptible fraction
I = 0.01  # Initial infectious fraction
R = 0.0  # Initial recovered fraction

# Time steps
num_steps = 1000
time_step = 0.1

# Lists to store simulation results
susceptible_list = [S]
infectious_list = [I]
recovered_list = [R]

# Simulation loop
for step in range(num_steps):
    # Calculate new fractions for each compartment
    new_S = S - beta * S * I * time_step + lambda_val * R * time_step
    new_I = I + (beta * S * I - gamma * I) * time_step
    new_R = R + gamma * I * time_step - lambda_val * R * time_step

    # Update compartment fractions
    S = new_S
    I = new_I
    R = new_R

    # Store results
    susceptible_list.append(S)
    infectious_list.append(I)
    recovered_list.append(R)

# Time points for plotting
t = [time_step * i for i in range(num_steps + 1)]

# Plot the results
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(t, susceptible_list, label='Susceptible')
plt.plot(t, infectious_list, label='Infectious')
plt.plot(t, recovered_list, label='Recovered')
plt.xlabel('Time')
plt.ylabel('Fraction of Population')
plt.title('SIRS Model Simulation (Homogeneously Mixed)')
plt.legend()
plt.grid()
plt.show()
