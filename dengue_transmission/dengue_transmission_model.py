from sys import maxsize

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.stats import alpha

# Definition of the gamma function (Rate of conversion from larvae to mosquitoes)
def gamma_t(t):
    day_of_year = t % 365
    if day_of_year < 120:  # Rainy season
        return 0.346
    elif day_of_year < 220:  # Normal season
        return 0.323
    else:  # Dry season
        return 0.091

# Adult mosquito mortality rate
def mu_1_t(t):
    day_of_year = t % 365
    if day_of_year < 120:  # Rainy season
        return 0.042
    elif day_of_year < 220:  # Normal season
        return 0.040
    else:  # Dry season
        return 0.059

# Mosquito population carrying capacity
def K0_t(t):
    day_of_year = t % 365
    if day_of_year < 120:  # Rainy season
        return 4500000
    elif day_of_year < 220:  # Normal season
        return 4000000
    else:  # Dry season
        return 3500000

# Mortality rate due to insecticide application
def mu_i_t(t):
    day_of_year = t % 365
    if day_of_year < 120:  # Rainy season
        return 0.958 * np.exp(-0.0151 * t)
    elif day_of_year < 220:  # Normal season
        return 0.960 * np.exp(-0.0151 * t)
    else:  # Dry season
        return 0.941 * np.exp(-0.0149 * t)

# Mortality rate due to larvicide application
def mu_l_t(t):
    day_of_year = t % 365
    if day_of_year < 120:  # Rainy season
        return 0.810 * np.exp(-0.00139 * t)
    elif day_of_year < 220:  # Normal season
        return 0.825 * np.exp(-0.00141 * t)
    else:  # Dry season
        return 0.884 * np.exp(-0.00145 * t)

# Model of the system of differential equations
def model_controls(y, t, params):
    A, MS, ME1, ME2, MI1, MI2, HS, HE1, HE2, HI1, HI2, HR1, HR2, HE12, HE21, HI12, HI21, HR = y
    # Unpack the parameters
    r, mu_2, sigma_M1, sigma_M2, mu_H, sigma_H1, sigma_H2, alpha_1, alpha_2, beta_M1, beta_M2, beta_H1, beta_H2, lambd, v, K_i, i, l = params
    # Seasonal parameters
    gamma = gamma_t(t)
    mu_1 = mu_1_t(t)
    K0 = K0_t(t)
    if i == 1:
        mu_i = mu_i_t(t)
    else: mu_i = 0
    if l == 1:
        mu_l = mu_l_t(t)
    else: mu_l = 0

    # Total populations
    M = MS + ME1 + ME2 + MI1 + MI2
    H = HS + HE1 + HE2 + HI1 + HI2 + HR1 + HR2 + HE12 + HE21 + HI12 + HI21 + HR
    # Differential equations
    dAdt = r * (1 - (A / (K_i * K0))) * M - (mu_2 + mu_l) * A - gamma * A
    dMSdT = gamma * A - (mu_1 + mu_i) * MS - beta_M1 * MS * ((HI1 + HI21) / H) - beta_M2 * MS * ((HI2 + HI12) / H)
    dME1dt = beta_M1 * MS * ((HI1 + HI21) / H) - (mu_1 + mu_i) * ME1 - sigma_M1 * ME1
    dME2dt = beta_M2 * MS * ((HI2 + HI12) / H) - (mu_2 + mu_i) * ME2 - sigma_M2 * ME2
    dMI1dt = sigma_M1 * ME1 - (mu_1 + mu_i) * MI1
    dMI2dt = sigma_M2 * ME2 - (mu_1 + mu_i) * MI2
    dHSdt = mu_H * H - mu_H * HS - beta_H1 * HS * (MI1 / M) - beta_H2 * HS * (MI2 / M) - v * HS
    dHE1dt = beta_H1 * HS * (MI1 / M) - sigma_H1 * HE1 - mu_H * HE1
    dHE2dt = beta_H2 * HS * (MI2 / M) - sigma_H2 * HE2 - mu_H * HE2
    dHI1dt = sigma_H1 * HE1 - mu_H * HI1 - alpha_1 * HI1
    dHI2dt = sigma_H2 * HE2 - mu_H * HI2 - alpha_2 * HI2
    dHR1dt = alpha_1 * HI1 - mu_H * HR1 - lambd * beta_H2 * HR1 * (MI2 / M)
    dHR2dt = alpha_2 * HI2 - mu_H * HR2 - lambd * beta_H1 * HR2 * (MI1 / M)
    dHE12dt = lambd * beta_H2 * HR1 * (MI2 / M) - mu_H * HE12 - sigma_H2 * HE12
    dHE21dt = lambd * beta_H1 * HR2 * (MI1 / M) - mu_H * HE21 - sigma_H1 * HE21
    dHI12dt = sigma_H2 * HE12 - mu_H * HI12 - alpha_2 * HI12
    dHI21dt = sigma_H1 * HE21 - mu_H * HI21 - alpha_1 * HI21
    dHRdt = alpha_2 * HI12 + alpha_1 * HI21 - mu_H * HR + v * HS
    return [dAdt, dMSdT, dME1dt, dME2dt, dMI1dt, dMI2dt, dHSdt, dHE1dt, dHE2dt, dHI1dt, dHI2dt, dHR1dt, dHR2dt, dHE12dt, dHE21dt, dHI12dt, dHI21dt, dHRdt]

def model_null(y, t, params):
    A, MS, ME1, ME2, MI1, MI2, HS, HE1, HE2, HI1, HI2, HR1, HR2, HE12, HE21, HI12, HI21, HR = y
    # Unpack the parameters
    r, mu_2, sigma_M1, sigma_M2, mu_H, sigma_H1, sigma_H2, alpha_1, alpha_2, beta_M1, beta_M2, beta_H1, beta_H2, lambd = params
    # Seasonal parameters
    gamma = gamma_t(t)
    mu_1 = mu_1_t(t)
    K0 = K0_t(t)

    # Total populations
    M = MS + ME1 + ME2 + MI1 + MI2
    H = HS + HE1 + HE2 + HI1 + HI2 + HR1 + HR2 + HE12 + HE21 + HI12 + HI21 + HR
    # Differential equations
    dAdt = r * (1 - (A / K0)) * M - mu_2 * A - gamma * A
    dMSdT = gamma * A - mu_1 * MS - beta_M1 * MS * ((HI1 + HI21) / H) - beta_M2 * MS * ((HI2 + HI12) / H)
    dME1dt = beta_M1 * MS * ((HI1 + HI21) / H) - mu_1 * ME1 - sigma_M1 * ME1
    dME2dt = beta_M2 * MS * ((HI2 + HI12) / H) - mu_2 * ME2 - sigma_M2 * ME2
    dMI1dt = sigma_M1 * ME1 - mu_1 * MI1
    dMI2dt = sigma_M2 * ME2 - mu_1 * MI2
    dHSdt = mu_H * H - mu_H * HS - beta_H1 * HS * (MI1 / M) - beta_H2 * HS * (MI2 / M)
    dHE1dt = beta_H1 * HS * (MI1 / M) - sigma_H1 * HE1 - mu_H * HE1
    dHE2dt = beta_H2 * HS * (MI2 / M) - sigma_H2 * HE2 - mu_H * HE2
    dHI1dt = sigma_H1 * HE1 - mu_H * HI1 - alpha_1 * HI1
    dHI2dt = sigma_H2 * HE2 - mu_H * HI2 - alpha_2 * HI2
    dHR1dt = alpha_1 * HI1 - mu_H * HR1 - lambd * beta_H2 * HR1 * (MI2 / M)
    dHR2dt = alpha_2 * HI2 - mu_H * HR2 - lambd * beta_H1 * HR2 * (MI1 / M)
    dHE12dt = lambd * beta_H2 * HR1 * (MI2 / M) - mu_H * HE12 - sigma_H2 * HE12
    dHE21dt = lambd * beta_H1 * HR2 * (MI1 / M) - mu_H * HE21 - sigma_H1 * HE21
    dHI12dt = sigma_H2 * HE12 - mu_H * HI12 - alpha_2 * HI12
    dHI21dt = sigma_H1 * HE21 - mu_H * HI21 - alpha_1 * HI21
    dHRdt = alpha_2 * HI12 + alpha_1 * HI21 - mu_H * HR
    return [dAdt, dMSdT, dME1dt, dME2dt, dMI1dt, dMI2dt, dHSdt, dHE1dt, dHE2dt, dHI1dt, dHI2dt, dHR1dt, dHR2dt, dHE12dt, dHE21dt, dHI12dt, dHI21dt, dHRdt]

# Parameters of the infected model
# r, mu_2, sigma_M1, sigma_M2, mu_H, sigma_H1, sigma_H2, alpha_1, alpha_2, beta_M1, beta_M2, beta_H1, beta_H2, lambd = params
params = [
    1,   # r: Intrinsic growth rate [1/t]
    0.05, # mu_2: Mortality rate of aquatic-stage mosquitoes [1/t]
    0.167,   # sigma_M1: Incubation period of serotype 1 in mosquitoes [1/t]
    0.25,   # sigma_M2: Incubation period of serotype 2 in mosquitoes [1/t]
    0.000042,   # mu_H: Human death and immigration rate [1/t]
    0.1,   # sigma_H1: Incubation period of serotype 1 in humans [1/t]
    0.142,  # sigma_H2: Incubation period of serotype 2 in humans [1/t]
    0.126,      # alpha_1: Recovery rate of humans from serotype 1 [1/t]
    0.142,      # alpha_2: Recovery rate of humans from serotype 2 [1/t]
    0.7,  # beta_M1: Effective contact rate between humans infected with serotype 1 and susceptible mosquitoes [1/t]
    0.75,   # beta_M2: Effective contact rate between humans infected with serotype 2 and susceptible mosquitoes [1/t]
    0.325,   # beta_H1: Effective contact rate between susceptible humans and mosquitoes infected with serotype 1 [1/t]
    0.375,    # beta_H2: Effective contact rate between susceptible humans and mosquitoes infected with serotype 2 [1/t]
    0.011    # lambda: Percentage reduction in secondary infection (cross-immunity) [U]
]

# Initial conditions
y0 = [2000000, 900000, 0, 0, 0, 0, 179668, 0, 0, 23, 55, 0, 0, 0, 0, 23, 55, 0]  # [A, MS, ME1, ME2, MI1, MI2, HS, HE1, HE2, HI1, HI2, HR1, HR2, HE12, HE21, HI12, HI21, HR]

# Simulation time
t = np.linspace(0, 365, 10000)

# Solve the system with null control
solution_null = odeint(model_null, y0, t, args=(params,))
A, MS, ME1, ME2, MI1, MI2, HS, HE1, HE2, HI1, HI2, HR1, HR2, HE12, HE21, HI12, HI21, HR = solution_null.T

# Solve the system with mechanical control
params_mech = params + [0, 0.5, 0, 0]
solution_mech = odeint(model_controls, y0, t, args=(params_mech,))
A_m, MS_m, ME1_m, ME2_m, MI1_m, MI2_m, HS_m, HE1_m, HE2_m, HI1_m, HI2_m, HR1_m, HR2_m, HE12_m, HE21_m, HI12_m, HI21_m, HR_m = solution_mech.T

# Solve the system with mechanical control and insecticide
params_mechins = params + [0, 0.5, 1, 0]
solution_mechins = odeint(model_controls, y0, t, args=(params_mechins,))
A_mi, MS_mi, ME1_mi, ME2_mi, MI1_mi, MI2_mi, HS_mi, HE1_mi, HE2_mi, HI1_mi, HI2_mi, HR1_mi, HR2_mi, HE12_mi, HE21_mi, HI12_mi, HI21_mi, HR_mi = solution_mechins.T

# Solve the system with mechanical control and larvicide
params_mechlar = params + [0, 0.5, 0, 1]
solution_mechlar = odeint(model_controls, y0, t, args=(params_mechlar,))
A_ml, MS_ml, ME1_ml, ME2_ml, MI1_ml, MI2_ml, HS_ml, HE1_ml, HE2_ml, HI1_ml, HI2_ml, HR1_ml, HR2_ml, HE12_ml, HE21_ml, HI12_ml, HI21_ml, HR_ml = solution_mechlar.T

# Solve the system with mechanical control and vaccination
params_mechvac = params + [0.013, 0.5, 0, 0]
solution_mechvac = odeint(model_controls, y0, t, args=(params_mechvac,))
A_mv, MS_mv, ME1_mv, ME2_mv, MI1_mv, MI2_mv, HS_mv, HE1_mv, HE2_mv, HI1_mv, HI2_mv, HR1_mv, HR2_mv, HE12_mv, HE21_mv, HI12_mv, HI21_mv, HR_mv = solution_mechvac.T

# Solve the system with insecticide and larvicide
params_inslar = params + [0, 1, 1, 1]
solution_inslar = odeint(model_controls, y0, t, args=(params_inslar,))
A_il, MS_il, ME1_il, ME2_il, MI1_il, MI2_il, HS_il, HE1_il, HE2_il, HI1_il, HI2_il, HR1_il, HR2_il, HE12_il, HE21_il, HI12_il, HI21_il, HR_il = solution_inslar.T

# Solve the system with insecticide and vaccination
params_insvac = params + [0.013, 1, 1, 0]
solution_insvac = odeint(model_controls, y0, t, args=(params_insvac,))
A_iv, MS_iv, ME1_iv, ME2_iv, MI1_iv, MI2_iv, HS_iv, HE1_iv, HE2_iv, HI1_iv, HI2_iv, HR1_iv, HR2_iv, HE12_iv, HE21_iv, HI12_iv, HI21_iv, HR_iv = solution_insvac.T

# Solve the system with larvicide
params_lar = params + [0, 1, 0, 1]
solution_lar = odeint(model_controls, y0, t, args=(params_lar,))
A_l, MS_l, ME1_l, ME2_l, MI1_l, MI2_l, HS_l, HE1_l, HE2_l, HI1_l, HI2_l, HR1_l, HR2_l, HE12_l, HE21_l, HI12_l, HI21_l, HR_l = solution_lar.T

# Solve the system with larvicide and vaccination
params_larvac = params + [0.013, 1, 0, 1]
solution_larvac = odeint(model_controls, y0, t, args=(params_larvac,))
A_lv, MS_lv, ME1_lv, ME2_lv, MI1_lv, MI2_lv, HS_lv, HE1_lv, HE2_lv, HI1_lv, HI2_lv, HR1_lv, HR2_lv, HE12_lv, HE21_lv, HI12_lv, HI21_lv, HR_lv = solution_larvac.T

# Solve the system with intensive vaccination
params_vac = params + [0.13, 1, 0, 0]
solution_vac = odeint(model_controls, y0, t, args=(params_vac,))
A_v, MS_v, ME1_v, ME2_v, MI1_v, MI2_v, HS_v, HE1_v, HE2_v, HI1_v, HI2_v, HR1_v, HR2_v, HE12_v, HE21_v, HI12_v, HI21_v, HR_v = solution_vac.T

# Plot results in separate graphs

# Graph 1: Healthy, Exposed, Infected, and Recovered humans with null control
plt.figure(figsize=(8, 5))
plt.plot(t, HS, label="HS")
plt.plot(t, HE1, label="HE1")
plt.plot(t, HE12, label="HE12")
plt.plot(t, HE2, label="HE2")
plt.plot(t, HE21, label="HE21")
plt.plot(t, HI1, label="HI1")
plt.plot(t, HI21, label="HI21")
plt.plot(t, HI2, label="HI2")
plt.plot(t, HI12, label="HI12")
plt.plot(t, HR1, label="HR1")
plt.plot(t, HR2, label="HR2")
plt.plot(t, HR, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with null control")
plt.grid()
plt.show()

# Graph 2: Healthy, Exposed, and Infected mosquitoes with null control
plt.figure(figsize=(8, 5))
plt.plot(t, MS, label="MS")
plt.plot(t, ME1, label="ME1")
plt.plot(t, ME2, label="ME2")
plt.plot(t, MI1, label="MI1")
plt.plot(t, MI2, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with null control")
plt.grid()
plt.show()

# Graph 3: Healthy, Exposed, Infected, and Recovered humans with mechanical control
plt.figure(figsize=(8, 5))
plt.plot(t, HS_m, label="HS")
plt.plot(t, HE1_m, label="HE1")
plt.plot(t, HE12_m, label="HE12")
plt.plot(t, HE2_m, label="HE2")
plt.plot(t, HE21_m, label="HE21")
plt.plot(t, HI1_m, label="HI1")
plt.plot(t, HI21_m, label="HI21")
plt.plot(t, HI2_m, label="HI2")
plt.plot(t, HI12_m, label="HI12")
plt.plot(t, HR1_m, label="HR1")
plt.plot(t, HR2_m, label="HR2")
plt.plot(t, HR_m, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with mechanical control")
plt.grid()
plt.show()

# Graph 4: Healthy, Exposed, and Infected mosquitoes with mechanical control
plt.figure(figsize=(8, 5))
plt.plot(t, MS_m, label="MS")
plt.plot(t, ME1_m, label="ME1")
plt.plot(t, ME2_m, label="ME2")
plt.plot(t, MI1_m, label="MI1")
plt.plot(t, MI2_m, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, 14500000, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, 14500000, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, 14500000, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with mechanical control")
plt.grid()
plt.show()

# Graph 5: Healthy, Exposed, Infected, and Recovered humans with mechanical control and insecticide
plt.figure(figsize=(8, 5))
plt.plot(t, HS_mi, label="HS")
plt.plot(t, HE1_mi, label="HE1")
plt.plot(t, HE12_mi, label="HE12")
plt.plot(t, HE2_mi, label="HE2")
plt.plot(t, HE21_mi, label="HE21")
plt.plot(t, HI1_mi, label="HI1")
plt.plot(t, HI21_mi, label="HI21")
plt.plot(t, HI2_mi, label="HI2")
plt.plot(t, HI12_mi, label="HI12")
plt.plot(t, HR1_mi, label="HR1")
plt.plot(t, HR2_mi, label="HR2")
plt.plot(t, HR_mi, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with mechanical control and insecticide")
plt.grid()
plt.show()

# Graph 6: Healthy, Exposed, and Infected mosquitoes with mechanical control and insecticide
plt.figure(figsize=(8, 5))
plt.plot(t, MS_mi, label="MS")
plt.plot(t, ME1_mi, label="ME1")
plt.plot(t, ME2_mi, label="ME2")
plt.plot(t, MI1_mi, label="MI1")
plt.plot(t, MI2_mi, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_mi) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_mi) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_mi) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with mechanical control and insecticide")
plt.grid()
plt.show()

# Graph 7: Healthy, Exposed, Infected, and Recovered humans with mechanical control and larvicide
plt.figure(figsize=(8, 5))
plt.plot(t, HS_ml, label="HS")
plt.plot(t, HE1_ml, label="HE1")
plt.plot(t, HE12_ml, label="HE12")
plt.plot(t, HE2_ml, label="HE2")
plt.plot(t, HE21_ml, label="HE21")
plt.plot(t, HI1_ml, label="HI1")
plt.plot(t, HI21_ml, label="HI21")
plt.plot(t, HI2_ml, label="HI2")
plt.plot(t, HI12_ml, label="HI12")
plt.plot(t, HR1_ml, label="HR1")
plt.plot(t, HR2_ml, label="HR2")
plt.plot(t, HR_ml, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_ml) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_ml) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_ml) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with mechanical control and larvicide")
plt.grid()
plt.show()

# Graph 8: Healthy, Exposed, and Infected mosquitoes with mechanical control and larvicide
plt.plot(t, MS_ml, label="MS")
plt.plot(t, ME1_ml, label="ME1")
plt.plot(t, ME2_ml, label="ME2")
plt.plot(t, MI1_ml, label="MI1")
plt.plot(t, MI2_ml, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_ml) * 0.90, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_ml) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_ml) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with mechanical control and larvicide")
plt.grid()
plt.show()

# Graph 9: Healthy, Exposed, Infected, and Recovered humans with mechanical control and vaccination
plt.figure(figsize=(8, 5))
plt.plot(t, HS_mv, label="HS")
plt.plot(t, HE1_mv, label="HE1")
plt.plot(t, HE12_mv, label="HE12")
plt.plot(t, HE2_mv, label="HE2")
plt.plot(t, HE21_mv, label="HE21")
plt.plot(t, HI1_mv, label="HI1")
plt.plot(t, HI21_mv, label="HI21")
plt.plot(t, HI2_mv, label="HI2")
plt.plot(t, HI12_mv, label="HI12")
plt.plot(t, HR1_mv, label="HR1")
plt.plot(t, HR2_mv, label="HR2")
plt.plot(t, HR_mv, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_mv) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_mv) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_mv) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with mechanical control and vaccination")
plt.grid()
plt.show()

# Graph 10: Healthy, Exposed, and Infected mosquitoes with mechanical control and vaccination
plt.plot(t, MS_mv, label="MS")
plt.plot(t, ME1_mv, label="ME1")
plt.plot(t, ME2_mv, label="ME2")
plt.plot(t, MI1_mv, label="MI1")
plt.plot(t, MI2_mv, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_mv) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_mv) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_mv) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with mechanical control and vaccination")
plt.grid()
plt.show()

# Graph 11: Healthy, Exposed, Infected, and Recovered humans with insecticide and larvicide
plt.figure(figsize=(8, 5))
plt.plot(t, HS_il, label="HS")
plt.plot(t, HE1_il, label="HE1")
plt.plot(t, HE12_il, label="HE12")
plt.plot(t, HE2_il, label="HE2")
plt.plot(t, HE21_il, label="HE21")
plt.plot(t, HI1_il, label="HI1")
plt.plot(t, HI21_il, label="HI21")
plt.plot(t, HI2_il, label="HI2")
plt.plot(t, HI12_il, label="HI12")
plt.plot(t, HR1_il, label="HR1")
plt.plot(t, HR2_il, label="HR2")
plt.plot(t, HR_il, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_il) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_il) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_il) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with insecticide and larvicide")
plt.grid()
plt.show()

# Graph 12: Healthy, Exposed, and Infected mosquitoes with insecticide and larvicide
plt.plot(t, MS_il, label="MS")
plt.plot(t, ME1_il, label="ME1")
plt.plot(t, ME2_il, label="ME2")
plt.plot(t, MI1_il, label="MI1")
plt.plot(t, MI2_il, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_il) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_il) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_il) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with insecticide and larvicide")
plt.grid()
plt.show()

# Graph 13: Healthy, Exposed, Infected, and Recovered humans with insecticide and vaccination
plt.figure(figsize=(8, 5))
plt.plot(t, HS_iv, label="HS")
plt.plot(t, HE1_iv, label="HE1")
plt.plot(t, HE12_iv, label="HE12")
plt.plot(t, HE2_iv, label="HE2")
plt.plot(t, HE21_iv, label="HE21")
plt.plot(t, HI1_iv, label="HI1")
plt.plot(t, HI21_iv, label="HI21")
plt.plot(t, HI2_iv, label="HI2")
plt.plot(t, HI12_iv, label="HI12")
plt.plot(t, HR1_iv, label="HR1")
plt.plot(t, HR2_iv, label="HR2")
plt.plot(t, HR_iv, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_iv) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_iv) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_iv) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with insecticide and vaccination")
plt.grid()
plt.show()

# Graph 14: Healthy, Exposed, and Infected mosquitoes with insecticide and vaccination
plt.plot(t, MS_iv, label="MS")
plt.plot(t, ME1_iv, label="ME1")
plt.plot(t, ME2_iv, label="ME2")
plt.plot(t, MI1_iv, label="MI1")
plt.plot(t, MI2_iv, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_iv) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_iv) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_iv) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with insecticide and vaccination")
plt.grid()
plt.show()

# Graph 15: Healthy, Exposed, Infected, and Recovered humans with larvicide
plt.figure(figsize=(8, 5))
plt.plot(t, HS_l, label="HS")
plt.plot(t, HE1_l, label="HE1")
plt.plot(t, HE12_l, label="HE12")
plt.plot(t, HE2_l, label="HE2")
plt.plot(t, HE21_l, label="HE21")
plt.plot(t, HI1_l, label="HI1")
plt.plot(t, HI21_l, label="HI21")
plt.plot(t, HI2_l, label="HI2")
plt.plot(t, HI12_l, label="HI12")
plt.plot(t, HR1_l, label="HR1")
plt.plot(t, HR2_l, label="HR2")
plt.plot(t, HR_l, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_l) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_l) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_l) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with larvicide")
plt.grid()
plt.show()

# Graph 16: Healthy, Exposed, and Infected mosquitoes with larvicide
plt.plot(t, MS_l, label="MS")
plt.plot(t, ME1_l, label="ME1")
plt.plot(t, ME2_l, label="ME2")
plt.plot(t, MI1_l, label="MI1")
plt.plot(t, MI2_l, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_l) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_l) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_l) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with larvicide and vaccination")
plt.grid()
plt.show()

# Graph 17: Healthy, Exposed, Infected, and Recovered humans with larvicide and vaccination
plt.figure(figsize=(8, 5))
plt.plot(t, HS_lv, label="HS")
plt.plot(t, HE1_lv, label="HE1")
plt.plot(t, HE12_lv, label="HE12")
plt.plot(t, HE2_lv, label="HE2")
plt.plot(t, HE21_lv, label="HE21")
plt.plot(t, HI1_lv, label="HI1")
plt.plot(t, HI21_lv, label="HI21")
plt.plot(t, HI2_lv, label="HI2")
plt.plot(t, HI12_lv, label="HI12")
plt.plot(t, HR1_lv, label="HR1")
plt.plot(t, HR2_lv, label="HR2")
plt.plot(t, HR_lv, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_lv) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_lv) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_lv) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with larvicide and vaccination")
plt.grid()
plt.show()

# Graph 17: Healthy, Exposed, and Infected mosquitoes with larvicide and vaccination
plt.plot(t, MS_lv, label="MS")
plt.plot(t, ME1_lv, label="ME1")
plt.plot(t, ME2_lv, label="ME2")
plt.plot(t, MI1_lv, label="MI1")
plt.plot(t, MI2_lv, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_lv) * 0.9, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_lv) * 0.9, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_lv) * 0.9, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with larvicide and vaccination")
plt.grid()
plt.show()

# Graph 18: Healthy, Exposed, Infected, and Recovered humans with intensive vaccination
plt.figure(figsize=(8, 5))
plt.plot(t, HS_v, label="HS")
plt.plot(t, HE1_v, label="HE1")
plt.plot(t, HE12_v, label="HE12")
plt.plot(t, HE2_v, label="HE2")
plt.plot(t, HE21_v, label="HE21")
plt.plot(t, HI1_v, label="HI1")
plt.plot(t, HI21_v, label="HI21")
plt.plot(t, HI2_v, label="HI2")
plt.plot(t, HI12_v, label="HI12")
plt.plot(t, HR1_v, label="HR1")
plt.plot(t, HR2_v, label="HR2")
plt.plot(t, HR_v, label="HR")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(HS_v) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(HS_v) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(HS_v) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model human population with intensive vaccination")
plt.grid()
plt.show()

# Graph 19: Healthy, Exposed, and Infected mosquitoes with intensive vaccination
plt.plot(t, MS_v, label="MS")
plt.plot(t, ME1_v, label="ME1")
plt.plot(t, ME2_v, label="ME2")
plt.plot(t, MI1_v, label="MI1")
plt.plot(t, MI2_v, label="MI2")
plt.axvline(120, color="red", linewidth=1, linestyle="dashed")
plt.axvline(220, color="red", linewidth=1, linestyle="dashed")
plt.text(60, max(MS_v) * 0.95, "Rainy", fontsize=12, bbox=dict(facecolor='lightblue', alpha=0.5))
plt.text(170, max(MS_v) * 0.95, "Normal", fontsize=12, bbox=dict(facecolor='lightgreen', alpha=0.5))
plt.text(270, max(MS_v) * 0.95, "Dry", fontsize=12, bbox=dict(facecolor='lightcoral', alpha=0.5))
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend()
plt.title("Dengue model mosquitoes with intensive vaccination")
plt.grid()
plt.show()
