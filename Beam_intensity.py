import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm
from scipy.integrate import simps
from scipy.special import wofz
import string


# Import data 
os.chdir('Beam intensity\FullRange')
data = pd.read_csv('crosssection.txt', sep='\t', decimal=',')  # Update file name/path as needed
E = data['Eg']  
sigma_data = data['sigma'] 
Dsigma = data['Dsigma']  

# Define functions
def compute_integral_with_uncertainty(sigma_fit, energy_distribution_values, E, uncertainty_factor=0.06):
    # Integral of the product with and without uncertainty corridor
    product_values = sigma_fit * energy_distribution_values
    sigma_integral = simps(product_values, E) * 10e-28  # Integral of the product
    
    # 6% uncertainty corridor
    sigma_upper = sigma_fit * (1 + uncertainty_factor)  # 6% higher
    sigma_lower = sigma_fit * (1 - uncertainty_factor)  # 6% lower
    
    product_upper = sigma_upper * energy_distribution_values
    product_lower = sigma_lower * energy_distribution_values
    
    integral_upper = simps(product_upper, E) * 10e-28
    integral_lower = simps(product_lower, E) * 10e-28
    
    uncertainty_integral = (integral_upper - integral_lower) / 2  # Uncertainty as half the difference
    
    return sigma_integral, uncertainty_integral
def voigt(E, sigma, gamma):
    z = ((E - E.mean()) + 1j * gamma) / (sigma * np.sqrt(2))
    z = np.array(z)  # Ensure z is a NumPy array
    return wofz(z).real / (sigma * np.sqrt(2 * np.pi))

# Global parameters
S_n = 6.512  
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


# Parameters for the skewed normal distribution
m = 9.85189
s = 10 * 27.52646648793018968603973917197436094284 / 1000
b = -1.8881705022515593572762782059726305306 / 1000

a = b / s  # Shape parameter for skewness
scale = s
loc = m

decay_constant = np.log(2) / half_life
decays = dN / (epsilon * abs_intensity)
atoms_per_area = N_A * x * rho / A  # atoms/cm^2
t_corr = (1 - np.exp(-decay_constant * t_IRR)) * np.exp(-decay_constant * t_delay) * (1 - np.exp(-decay_constant * t_measurement))

# Uncertainties
u_half_life = 0.0006 * 24 * 3600  # seconds
u_r = 0.01 * 2.54  # uncertainty in radius
u_S = u_r * np.pi * 1.05 / 2  # uncertainty in area
u_m = 1 * 10 ** (-3)  # uncertainty in mass
u_dN = np.array([72, 132, 37])  # uncertainties in dN
u_epsilon = np.array([0.0011, 0.0011, 0.0009])  # uncertainties in epsilon
u_abs_intensity = np.array([0.0095, 0.03, 0.003])  # uncertainties in intensity

u_decays = np.sqrt((u_dN / (epsilon * abs_intensity)) ** 2 +
                   (u_epsilon * dN / (epsilon ** 2 * abs_intensity)) ** 2 +
                   (u_abs_intensity * dN / (epsilon * abs_intensity ** 2)) ** 2)    


# Define the models function for fitting

# 1. Polynomial with increased number of terms
def model_polynomial_increased(E, sigma_0, p, a0, a1, a2, a3, a4):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    poly = a0 + a1 * E + a2 * E**2 + a3 * E**3 + a4 * E**4
    return base * poly

# 2. Polynomial with decreased number of terms
def model_polynomial_decreased(E, sigma_0, p, a0, a1):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    poly = a0 + a1 * E
    return base * poly

# 3. Voigt profile
def model_voigt(E, sigma_0, p, sigma, gamma):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base * voigt(E.to_numpy(), sigma, gamma)  # Ensure E is a NumPy array

# 4. Increased power-law
def model_power_increased(E, sigma_0, p, a):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base * E**a

# 5. Decreased power-law
def model_power_decreased(E, sigma_0, p, a):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base / E**a

# 6. Sinh function
def model_sinh(E, sigma_0, p, A):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base * np.sinh(A * E)

# 7. Cosh function
def model_cosh(E, sigma_0, p, A):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base * np.cosh(A * E)

# 8. Original model
def original_model(E, sigma_0, p, a0, a1, a2, a3):
    base = sigma_0 * ((E - S_n) / S_n) ** p
    return base * (a0 + a1 * E + a2 * E ** 2 + a3 * E ** 3)

# Define all models in a dictionary
models = {
    "Polynomial (Increased)": model_polynomial_increased,
    "Polynomial (Decreased)": model_polynomial_decreased,
    "Voigt Profile": model_voigt,
    "Power Law (Increased)": model_power_increased,
    "Power Law (Decreased)": model_power_decreased,
    "Sinh": model_sinh,
    "Cosh": model_cosh,
    "Original Model": original_model
}

# Initial guesses for each model
initial_guesses = {
    "Polynomial (Increased)": [1, 1, 1, 1, 1, 1, 1],
    "Polynomial (Decreased)": [1, 1, 1, 1],
    "Voigt Profile": [1, 1, 1, 1],
    "Power Law (Increased)": [1, 1, 1],
    "Power Law (Decreased)": [1, 1, 1],
    "Sinh": [1, 1, 1],
    "Cosh": [1, 1, 1],
    "Original Model": [1, 1, 1, 1, 1, 1]

}
latex_results = []

# Fit each model and plot results
for model_name, model_func in models.items():
    try:
        # Fit the model
        popt, pcov = curve_fit(
            model_func,
            E,
            sigma_data,
            p0=initial_guesses[model_name],
            sigma=Dsigma,
            maxfev=10000,  # Increase iterations if needed
        )
        
        # Calculate fitted values and residuals
        sigma_fit = model_func(E, *popt)
        residuals = sigma_data - sigma_fit
        
        # Print parameters
        print(f"\n{model_name} - Fitted Parameters:")
        for i, param in enumerate(popt):
            print(f"Param {i + 1}: {param:.4f}")
            
        chi_squared = np.sum((residuals ** 2) / Dsigma ** 2)
        
        fitted_parameters = popt
        N = len(sigma_data)  # Number of data points
        p = len(fitted_parameters)  # Number of fitted parameters
        chi_squared_reduced = chi_squared / (N - p)

        print(f"\nChi-squared: χ² = {chi_squared:.4f}")
        print(f"Reduced Chi-squared: χ²ν = {chi_squared_reduced:.4f}")
        
        # Calculate the uncertainty corridor (6%)
        uncertainty_factor = 0.06
        sigma_upper = sigma_fit * (1 + uncertainty_factor)  # 6% higher
        sigma_lower = sigma_fit * (1 - uncertainty_factor)  # 6% lower
        
        energy_distribution_values = skewnorm.pdf(E, a, loc=loc, scale=scale)
        sigma_integral, sigma_uncertainty_integral = computze_integral_with_uncertainty(sigma_fit, energy_distribution_values, E, uncertainty_factor)
        
        
        
        I_sigma = decays * decay_constant / (t_corr * atoms_per_area)
        u_I_sigma = np.sqrt((u_decays * decay_constant / (t_corr * atoms_per_area)) ** 2)

    
        for i in range(len(I_sigma)):
            print(f"I_sigma[{i}] = {I_sigma[i]:.3g} ± {u_I_sigma[i]:.2g} (1/s)")

        # Weighted average of intensity
        weights = 1 / u_I_sigma ** 2
        weighted_I_sigma = np.sum(I_sigma * weights) / np.sum(weights)
        u_weighted_I_sigma = 1 / np.sqrt(np.sum(weights))

        print(f"Mean: I_sigma = {weighted_I_sigma:.3g} ± {u_weighted_I_sigma:.2g} (barn/s)")
        
        I=weighted_I_sigma/sigma_integral
        u_I= sigma_uncertainty_integral*weighted_I_sigma/sigma_integral**2
        u_I_relative=100*u_I/I
        # Printing final results with updated formatting
        print(f"Final result: I = {I:.2g} ± {u_I:.1g} (Hz), relative uncertainty: {u_I_relative:.2g}%")


        # Fitted parameters and chi-squared values for the legend
        param_str = ', '.join([f'{string.ascii_lowercase[i]}: {param:.4f}' for i, param in enumerate(popt)])
        chi_squared_str = f"Chi-squared: χ² = {chi_squared:.4f}"
        chi_squared_reduced_str = f"Reduced χ²ν = {chi_squared_reduced:.4f}"
        
        # Prepare the legend string with fitted parameters and chi-squared values
        legend_str = f'{model_name} - Fitted Curve\n{param_str}\n{chi_squared_str}\n{chi_squared_reduced_str}'
        
        # Plot the fitted curve with updated legend
        plt.figure(figsize=(10, 6))
        plt.errorbar(E, sigma_data, yerr=Dsigma, fmt='o', label='Data', color='blue')

        plt.fill_between(E, sigma_lower, sigma_upper, color='gray', alpha=0.3, label='Uncertainty Corridor (±6%)')
        plt.plot(E, sigma_fit, 'r-', label=legend_str)  # Add the legend string here
        plt.title(f'{model_name} - Fitted Curve')
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Cross-Section (mbarn)')
        plt.legend()
        plt.grid(True)


   
        # Save the plot as PNG
        plt.savefig(f'{model_name}_fitted_curve.png', dpi=1000)
        plt.show()
        
        # Plot residuals
        plt.figure(figsize=(10, 6))
        plt.errorbar(E, residuals, yerr=Dsigma, fmt='o', label='Residuals', color='orange')
        plt.axhline(0, color='red', linestyle='--', label='Zero Line')
        plt.title(f'{model_name} - Residual Plot')
        plt.xlabel('Energy (MeV)')
        plt.ylabel('Residuals (mbarn)')
        plt.legend()
        plt.grid(True)
        
        # Save the residual plot as PNG
        plt.savefig(f'{model_name}_residual_plot.png', dpi=1000)
        plt.show()
        
        latex_results.append({
            "model": model_name,
            "chi_squared": chi_squared,
            "chi_squared_reduced": chi_squared_reduced,
            "sigma": sigma_integral,
            "sigma_uncertainty": sigma_uncertainty_integral,
            "I": I,
            "I_uncertainty": u_I,
            "Relative uncertainty": u_I_relative
        })


    
    except Exception as e:
        print(f"\n{model_name} - Error during fitting: {e}")

# Adjust LaTeX table formatting for significant digits
latex_table = r"\begin{table}[h!]\centering" + "\n"
latex_table += r"\begin{tabular}{|l|c|c|c|c|c|c|c|}" + "\n"
latex_table += r"\hline\nModel & $\chi^2$ & $\chi^2_\nu$ & $\sigma$ & $\Delta\sigma$ & $I$ & $\Delta I$ & $\Delta I / I$  \\\hline" + "\n"

for result in latex_results:
    row_color = r"\rowcolor{green!20}" if result["chi_squared_reduced"] < 1 else ""
    latex_table += (
        f"{row_color} {result['model']} & "
        f"{result['chi_squared']:.2g} & "
        f"{result['chi_squared_reduced']:.2g} & "
        f"{result['sigma']:.3g} & "
        f"{result['sigma_uncertainty']:.2g} & "
        f"{result['I']:.3g} & "
        f"{result['I_uncertainty']:.2g} & "
        f"{result['Relative uncertainty']:.2g}\\% \\\\\n"
    )

latex_table += r"\hline\n\end{tabular}" + "\n"
latex_table += r"\caption{Fit Results with $\chi^2$, Reduced $\chi^2$, and Relative Uncertainty. Rows with $\chi^2_\nu < 1$ are highlighted.}" + "\n"
latex_table += r"\end{table}"

# Save LaTeX table to a file
latex_file_path = "fit_results_table.tex"
with open(latex_file_path, "w") as f:
    f.write(latex_table)

print("LaTeX table saved to 'fit_results_table.tex'.")
