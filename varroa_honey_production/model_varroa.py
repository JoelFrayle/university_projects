import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Differential equations system model
def model(y, t, params):
    # Unpack the variables
    C, PH, PM, B, BI, A, I, M = y
    # Unpack the parameters
    gamma, mu, phi, theta, tau, omega, rho, beta, alpha, K, delta, lamda, lamda2 = params

    # Differential equations
    dCdt = beta * PH - (I / (I + A)) * C  # Dynamics of healthy cells
    dPHdt = (I / (I + A)) * C - beta * PH  # Dynamics of male parasites
    dPMdt = alpha * PH * delta - beta * PM  # Dynamics of female parasites
    dBdt = (gamma * (1 - ((A + I) / K)) - gamma * (PH / (PH + C)) * (1 - ((A + I) / K))
            - delta * B - lamda * B)  # Dynamics of healthy larvae
    dBIdt = (gamma * (PH / (PH + C)) * (1 - ((A + I) / K)) - delta * BI - lamda2 * BI)  # Dynamics of infected larvae
    dAdt = delta * B - phi * I - mu * A  # Dynamics of healthy adult bees
    dIdt = delta * BI + phi * I - omega * I  # Dynamics of infected adult bees
    dMdt = theta * A + tau * I  # Honey production

    return [dCdt, dPHdt, dPMdt, dBdt, dBIdt, dAdt, dIdt, dMdt]


# Parameters for the infected colony model
params = [
    2500,   # gamma: Daily birth rate of bees
    0.04,   # mu: Mortality rate of healthy bees
    0.04959,  # phi: Infection rate
    0.4,    # theta: Honey production rate by healthy bees
    0.02,   # tau: Honey production rate by infected bees
    0.2,    # omega: Mortality rate of infected bees
    0.02,   # rho: Birth rate of parasites
    0.007,  # beta: Mortality rate of parasites
    6,      # alpha: Female parasites per male
    80000,  # K: Carrying capacity
    1/21,   # delta: Maturation rate of larvae
    0.02,   # lamda: Mortality rate of healthy larvae
    0.05    # lamda2: Mortality rate of infected larvae
]

# Initial conditions
y0_infected = [500, 50, 30, 2500, 10, 50000, 10, 0]  # Initial values for the infected colony
y0_healthy = [500, 0, 0, 2500, 0, 50000, 0, 0]       # Initial values for the healthy colony

# Simulation time
t = np.linspace(0, 100, 1000)  # 100 days, 1000 time points

# Solve the system for the infected colony
solution_infected = odeint(model, y0_infected, t, args=(params,))
C, PH, PM, B, BI, A, I, M = solution_infected.T

# Solve the system for the healthy colony
solution_healthy = odeint(model, y0_healthy, t, args=(params,))
C_healthy, PH_healthy, PM_healthy, B_healthy, BI_healthy, A_healthy, I_healthy, M_healthy = solution_healthy.T

# Plot results in separate graphs

# Plot 1: Healthy and infected larvae
plt.figure(figsize=(8, 5))
plt.plot(t, B, label="Healthy larvae (B)", linestyle="--")
plt.plot(t, BI, label="Infected larvae (BI)", linestyle=":")
plt.xlabel("Time (days)")
plt.ylabel("Population")
plt.legend()
plt.title("Healthy and infected larvae")
plt.grid()
plt.show()

# Plot 2: Male and female parasites
plt.figure(figsize=(8, 5))
plt.plot(t, PH, label="Male parasites (PH)", linestyle="--")
plt.plot(t, PM, label="Female parasites (PM)", linestyle=":")
plt.xlabel("Time (days)")
plt.ylabel("Population")
plt.legend()
plt.title("Male and female parasites")
plt.grid()
plt.show()

# Plot 3: Healthy and infected adult bees
plt.figure(figsize=(8, 5))
plt.plot(t, A, label="Healthy bees in infected colony (A)", linestyle="--")
plt.plot(t, I, label="Infected bees (I)", linestyle=":")
plt.plot(t, A_healthy, label="Bees in healthy colony", linestyle="-")
plt.xlabel("Time (days)")
plt.ylabel("Population")
plt.legend()
plt.title("Adult bees (healthy and infected)")
plt.grid()
plt.show()

# Plot 4: Honey production
plt.figure(figsize=(8, 5))
plt.plot(t, M, label="Honey (infected colony)", linestyle="--")
plt.plot(t, M_healthy, label="Honey (healthy colony)", linestyle=":")
plt.xlabel("Time (days)")
plt.ylabel("Honey production (mL)")
plt.legend()
plt.title("Honey production in healthy and infected colonies")
plt.grid()
plt.show()
