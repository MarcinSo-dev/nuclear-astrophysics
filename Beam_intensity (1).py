import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm
from scipy.integrate import simps


# Cross-section fitting
data = pd.read_csv('crosssection.txt', sep='\t', decimal=',')  # Update file name/path as needed
E = data['Eg']  # Energy in keV
sigma_data = data['sigma']  # Cross-section in mbarn
Dsigma = data['Dsigma']  # Uncertainties in cross-section

# Define the model function for fitting
def cross_section_model(E, sigma_0, p, a0, a1, a2, a3):
    S_n = 6.512  # Given constant
    return sigma_0 * ((E - S_n) / S_n) ** p * (a0 + a1 * E + a2 * E ** 2 + a3 * E ** 3)

# Initial guess for the parameters [sigma_0, p, a0, a1, a2, a3]
initial_guess = [1, 1, 1, 1, 1, 1]

# Fit the model to the data
popt, pcov = curve_fit(cross_section_model, E, sigma_data, p0=initial_guess, sigma=Dsigma)

# Extract the fitted parameters and their errors
fitted_parameters = popt
parameter_errors = np.sqrt(np.diag(pcov))

# Print fitted parameters
param_names = ['sigma_0', 'p', 'a0', 'a1', 'a2', 'a3']
print("\nFitted parameters:")
for name, value, error in zip(param_names, fitted_parameters, parameter_errors):
    print(f"{name} = {value:.4f} ± {error:.4f}")

# Calculate fitted values and residuals
sigma_fit = cross_section_model(E, *popt)
residuals = sigma_data - sigma_fit

# Calculate chi-squared and reduced chi-squared
chi_squared = np.sum((residuals ** 2) / Dsigma ** 2)
N = len(sigma_data)  # Number of data points
p = len(fitted_parameters)  # Number of fitted parameters
chi_squared_reduced = chi_squared / (N - p)

print(f"\nChi-squared: χ² = {chi_squared:.4f}")
print(f"Reduced Chi-squared: χ²ν = {chi_squared_reduced:.4f}")

# Parameters for the skewed normal distribution
m = 9.85189
s = 10 * 27.52646648793018968603973917197436094284 / 1000
b = -1.8881705022515593572762782059726305306 / 1000

a = b / s  # Shape parameter for skewness
scale = s
loc = m

# Generate x values for the skewed normal distribution
x_dist = np.linspace(m - 4 * s, m + 4 * s, 10000)

# Generate skewed normal distribution values
pdf = skewnorm.pdf(x_dist, a, loc=loc, scale=scale)

# Calculate the energy distribution values for the energies in the cross-section data
energy_distribution_values = skewnorm.pdf(E, a, loc=loc, scale=scale)

# Step 3: Calculate the average product of the cross-section and energy distribution
product_values = sigma_fit * energy_distribution_values
sigma = simps(product_values, E)*10e-28 # Integral of the product

# Display the results
print(f'Cross section: {sigma} mbarn')

# Plot the original data and the fitted curve
plt.figure(figsize=(10, 6))
plt.errorbar(E, sigma_data, yerr=Dsigma, fmt='o', label='Data', color='blue')
plt.plot(E, sigma_fit, 'r-', label='Fitted Curve')
plt.title('Cross-Section Data and Fitted Curve')
plt.xlabel('Energy (keV)')
plt.ylabel('Cross-Section (mbarn)')
plt.legend()
plt.grid(True)

plt.show()

# Residual Plot
plt.figure(figsize=(10, 6))
plt.errorbar(E, residuals, yerr=Dsigma, fmt='o', label='Residuals', color='orange')
plt.axhline(0, color='red', linestyle='--', label='Zero Line')
plt.title('Residual Plot')
plt.xlabel('Energy (keV)')
plt.ylabel('Residuals (mbarn)')
plt.legend()
plt.grid(True)
plt.show()

# Plot the skewed normal distribution
plt.figure(figsize=(10, 6))
plt.plot(x_dist, pdf*2*10e23, label='Skewed Normal Distribution')  # Scaled for better visualization
plt.errorbar(E, sigma_data * 10 ** 22, yerr=Dsigma, fmt='o', label='Data', color='blue')

plt.title('Skewed Normal Distribution with Data')
plt.xlabel('x')
plt.ylabel('Probability Density')

plt.legend()
plt.grid(True)

plt.show()


# Constants & uncertainties
N_A = 6.023 * 10 ** 23  # atoms per mole
A = 197  # nuclei
rho = 19.283  # g/cm^3
x = 25.4 * 10 ** (-4)  # cm
t_IRR = 59 * 60  # seconds
t_delay = 19 * 60  # seconds
t_measurement = 12 * 3600  # seconds

half_life = 6.1669 * 24 * 3600  # seconds
S = np.pi * (0.75 * 2.54 / 2) ** 2  # cm^2
m = 147 * 10 ** (-3)  # g
dN = np.array([4636, 16893, 1089])
epsilon = np.array([0.0228, 0.0216, 0.0187])
abs_intensity = np.array([0.2288, 0.87, 0.066])

# Uncertainties
u_half_life = 0.0006 * 24 * 3600  # seconds
u_r = 0.01 * 2.54  # uncertainty in radius
u_S = u_r * np.pi * 1.05 / 2  # uncertainty in area
u_m = 1 * 10 ** (-3)  # uncertainty in mass
u_dN = np.array([72, 132, 37])  # uncertainties in dN
u_epsilon = np.array([0.0011, 0.0011, 0.0009])  # uncertainties in epsilon
u_abs_intensity = np.array([0.0095, 0.03, 0.003])  # uncertainties in intensity

# Calculating decay constant and decays
decay_constant = np.log(2) / half_life
decays = dN / (epsilon * abs_intensity)
atoms_per_area = N_A * x * rho / A  # atoms/cm^2
t_corr = (1 - np.exp(-decay_constant * t_IRR)) * np.exp(-decay_constant * t_delay) * (1 - np.exp(-decay_constant * t_measurement))

# Uncertainties in decays
u_dc = u_half_life * np.log(2) / half_life ** 2
u_decays = np.sqrt((u_dN / (epsilon * abs_intensity)) ** 2 +
                   (u_epsilon * dN / (epsilon ** 2 * abs_intensity)) ** 2 +
                   (u_abs_intensity * dN / (epsilon * abs_intensity ** 2)) ** 2)

# Calculate intensity
I = decays * decay_constant / (sigma * t_corr * atoms_per_area)
u_I = np.sqrt((u_decays * decay_constant / (sigma * t_corr * atoms_per_area)) ** 2)

for i in range(len(I)):
    print(f"I[{i}] = {I[i]:.3g} ± {u_I[i]:.2g} (1/s)")

# Weighted average of intensity
weights = 1 / u_I ** 2
weighted_I = np.sum(I * weights) / np.sum(weights)
u_weighted_I = 1 / np.sqrt(np.sum(weights))

print(f"Final result: I = {weighted_I:.3g} ± {u_weighted_I:.2g} (1/s)")
