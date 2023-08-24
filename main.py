# Parameters
beta = 0.3  # Constant transmission rate
gamma = 0.1  # Recovery rate
lambda_val = 0.05  # Immunity loss rate
alpha = 0.08  # Immunity waning rate
sigma = 0.005  # Vaccination rate (fraction of susceptibles vaccinated per day)
initial_population = 1000  # Total population size
T_inf = 1.0 / gamma  # Mean infectious duration
T_imm = 1.0 / alpha  # Mean immunity duration

# Initial conditions
S = 0.99 * initial_population  # Initial susceptible count
I = 0.01 * initial_population  # Initial infectious count
R = 0.0  # Initial recovered count
V = 0.0  # Initial vaccinated count

# Time steps
num_steps = 1000
time_step = 0.1

# Lists to store simulation results
susceptible_list = [S]
infectious_list = [I]
recovered_list = [R]
vaccinated_list = [V]

# Simulation loop
for step in range(num_steps):
    # Calculate new counts for each compartment
    new_S = S - (beta * S * I / initial_population) * time_step + lambda_val * R * time_step - sigma * S * time_step
    new_I = I + ((beta * S * I / initial_population) - gamma * I) * time_step - I * (time_step / T_inf)
    new_R = R + (gamma * I - lambda_val * R) * time_step + I * (time_step / T_inf) - R * (time_step / T_imm)
    new_V = V + sigma * S * time_step + R * (time_step / T_imm)

    # Update compartment counts
    S = new_S
    I = new_I
    R = new_R
    V = new_V

    # Store results
    susceptible_list.append(S)
    infectious_list.append(I)
    recovered_list.append(R)
    vaccinated_list.append(V)

# Time points for plotting
t = [time_step * i for i in range(num_steps + 1)]

# Plot the results
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(t, susceptible_list, label='Susceptible')
plt.plot(t, infectious_list, label='Infectious')
plt.plot(t, recovered_list, label='Recovered')
plt.plot(t, vaccinated_list, label='Vaccinated')
plt.xlabel('Time')
plt.ylabel('Count')
plt.title('SIRS Model Simulation with Vaccination and Exponential Dynamics')
plt.legend()
plt.grid()
plt.show()
