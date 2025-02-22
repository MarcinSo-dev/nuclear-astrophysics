import numpy as np
from scipy.stats import skewnorm, cauchy
import matplotlib.pyplot as plt
from scipy.integrate import simpson as simps

def generate_normalized_lorentzian_data(E, E0, Gamma, beam_energy, Qvalue, num_samples=10000):
    # Generate random data from the Lorentzian distribution using the Cauchy distribution
    random_data = cauchy.rvs(loc=E0, scale=Gamma, size=num_samples)
    
    # Now, we normalize the random data to make sure the total area under the distribution is 1
    hist, bin_edges = np.histogram(random_data, bins=E, density=True)

    # Normalize the histogram (so the integral of the PDF is 1)
    hist /= simps(hist, dx=np.diff(bin_edges)[0])  # Use Simpson's rule to normalize the histogram
    
    # Interpolate beam energy to match the bin edges (to ensure compatibility)
    beam_energy_interp = np.interp(bin_edges[:-1], E, beam_energy)
    
    # Multiply by beam energy (interpolated to match the bin edges)
    hist = hist * beam_energy_interp
    
    return hist, bin_edges+Qvalue

def g(E, E0, Gamma, beam_energy, Qvalue, num_samples=10000):
    # Generate random data from the Lorentzian distribution using the Cauchy distribution
    random_data = cauchy.rvs(loc=E0, scale=Gamma, size=num_samples)
    
    # Now, we normalize the random data to make sure the total area under the distribution is 1
    hist, bin_edges = np.histogram(random_data, bins=E, density=True)

    # Normalize the histogram (so the integral of the PDF is 1)
    hist /= simps(hist, dx=np.diff(bin_edges)[0])  # Use Simpson's rule to normalize the histogram
    
    # Interpolate beam energy to match the bin edges (to ensure compatibility)
    beam_energy_interp = np.interp(bin_edges[:-1], E, beam_energy)
    
    # Multiply by beam energy (interpolated to match the bin edges)
    hist = hist
    return hist, bin_edges

dalton_to_MeV = 931.49410372
E = np.linspace(0, 15, 10000)

m = 9.85189
s = 10 * 27.52646648793018968603973917197436094284 / 1000
b = -1.8881705022515593572762782059726305306 / 1000

a = b / s  # Shape parameter for skewness
scale = s
loc = m

beam_energy = skewnorm.pdf(E, a, loc=loc, scale=scale)

O_16_mass_excess = -4.73700223 # MeV
C_12_mass_excess = 0 # MeV
Be_8_mass_excess = 4.941674 # MeV
alpha_mass_excess = 2.4249158715 # MeV

O_16_mass = 16*dalton_to_MeV+O_16_mass_excess
C_12_mass = 12*dalton_to_MeV+C_12_mass_excess
Be_8_mass = 8*dalton_to_MeV+Be_8_mass_excess
alpha_mass = 4*dalton_to_MeV+alpha_mass_excess

O_16_Qvalue = O_16_mass_excess - (C_12_mass_excess+alpha_mass_excess)
C_12_Qvalue = C_12_mass_excess - (Be_8_mass_excess+alpha_mass_excess)
Be_8_Qvalue = Be_8_mass_excess - (2*alpha_mass_excess)

O_resonance_energy =  9.84455 #MeV
O_resonance_gamma = 9.88 *10**(-9) #MeV

C_resonance_energy = 9.870 #MeV
C_resonance_gamma = 0.850

Be_resonance_energy = 3.031 #MeV
Be_resonance_gamma = 1.52315 #MeV

Be_ground_state_energy = 0 #MeV
Be_ground_state_gamma = 5.57 *10**(-3) #MeV

# Generate normalized Lorentzian distributions for each resonance, multiplied by beam energy
O_hist, O_bin_edges = generate_normalized_lorentzian_data(E, O_resonance_energy, O_resonance_gamma, beam_energy, O_16_Qvalue)
C_hist, C_bin_edges = generate_normalized_lorentzian_data(E, C_resonance_energy, C_resonance_gamma, beam_energy, C_12_Qvalue)
Be_hist, Be_bin_edges = g(E, Be_resonance_energy, Be_resonance_gamma, beam_energy, Be_8_Qvalue)
Be_gs_hist, Be_gs_bin_edges = g(E, Be_ground_state_energy, Be_ground_state_gamma, beam_energy, Be_8_Qvalue)

O_hist /= np.sum(O_hist)
C_hist /= np.sum(C_hist)
# Plot the normalized distributions multiplied by beam energy
plt.plot(O_bin_edges[:-1], O_hist, label="O Resonance")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('Energy in CMS distribution (16O)')
plt.xlim(2.6, 2.8)
plt.savefig('plot', dpi=500)
plt.show()

plt.plot(C_bin_edges[:-1], C_hist, label="C Resonance")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('Energy in CMS distribution (12C)')
plt.xlim(2, 2.6)
plt.show()

O_reaction_alpha_mass = O_bin_edges[:-1]*C_12_mass/(C_12_mass+alpha_mass)
O_reaction_C_mass = O_bin_edges[:-1]*alpha_mass/(C_12_mass+alpha_mass)

plt.plot(O_reaction_alpha_mass, O_hist, label = "16O reaction alpha energy")
plt.plot(O_reaction_C_mass, O_hist, label = "16O reaction 12C energy")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('O16 reaction nuclei energies')
plt.xlim(0, 2.5)
plt.show()


C_reaction_Be_ECM = C_bin_edges[:-1] * Be_8_mass / (C_12_mass + Be_8_mass)
C_hist_multiplication_ECM = C_hist * C_reaction_Be_ECM

# Plot the result
plt.plot(C_reaction_Be_ECM, C_hist_multiplication_ECM, label="C12  8Be Resonance ECM")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Multiplied Probability Density')
plt.title('C12 x Be Resonance Distribution')
plt.xlim(0, 2.5)
plt.show()

plt.plot(Be_bin_edges[:-1], Be_hist/10, label = "8Be_resonance")
plt.plot(C_bin_edges[:-1], C_hist, label="C Resonance")
plt.xlim(1, 3)
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('Energy in CMS distribution (12C)')
plt.show()


# Multiply the C and Be histograms element-wise
multiplied_hist = C_hist * Be_hist

# Plot the result
plt.plot(C_bin_edges[:-1], multiplied_hist, label="C Resonance x Be Resonance")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Multiplied Probability Density')
plt.title('Product of C12 and 8Be Resonance Distributions')
plt.xlim(1, 2.5)
plt.show()


# Multiply the C and Be histograms element-wise
multiplied_hist = C_hist * Be_gs_hist

# Plot the result
plt.plot(C_bin_edges[:-1], multiplied_hist, label="C Resonance x Be Resonance")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Multiplied Probability Density')
plt.title('Product of C12 and 8Be Resonance Distributions')
plt.xlim(0, 2.5)
plt.show()




def lorentzian(E, E0, Gamma):
    return (1 / np.pi) * (Gamma / ((E - E0)**2 + Gamma**2))

# Calculate the Lorentzian for Be resonance and ground state (without beam energy)
Be_resonance_lorentzian = lorentzian(E, Be_resonance_energy, Be_resonance_gamma)
Be_ground_state_lorentzian = lorentzian(E, Be_ground_state_energy, Be_ground_state_gamma)
# Trim the longer Lorentzian distributions to match the length of C_hist
Be_resonance_lorentzian = Be_resonance_lorentzian[:len(C_hist)]
Be_ground_state_lorentzian = Be_ground_state_lorentzian[:len(C_hist)]

# Adjust C_hist based on the resonance and ground state contributions
C_hist_adjusted = C_hist * (Be_resonance_lorentzian + Be_ground_state_lorentzian)

# Multiply the adjusted C_hist with Be_hist (which already includes beam energy considerations)
multiplied_hist_adjusted = C_hist_adjusted * Be_hist

# Plot the result
plt.plot(C_bin_edges[:-1], multiplied_hist_adjusted, label="Adjusted C Resonance x Be Resonance")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Multiplied Probability Density')
plt.title('Product of Adjusted C12 and Be8 Resonance Distributions (without beam energy)')
plt.xlim(2, 2.5)
plt.show()


C_reaction_alpha_mass = C_bin_edges[:-1] * Be_8_mass / (Be_8_mass + alpha_mass)
C_reaction_Be_mass = C_bin_edges[:-1] * alpha_mass / (Be_8_mass + alpha_mass)

plt.plot(C_reaction_alpha_mass, multiplied_hist_adjusted, label = "12C reaction alpha energy")
plt.plot(C_reaction_Be_mass, multiplied_hist_adjusted, label = "12C reaction 8Be energy")
plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('C12 reaction nuclei energies')
plt.xlim(0, 2.5)
plt.show()













from scipy.stats import norm

# Fit Gaussian to the C_reaction_alpha_mass and C_reaction_Be_mass data
mu_alpha, sigma_alpha = norm.fit(C_reaction_alpha_mass)
mu_Be, sigma_Be = norm.fit(C_reaction_Be_mass)

# Print the fitted parameters for both Gaussian distributions
print(f"Gaussian fit for C_reaction_alpha_mass:")
print(f"Mean: {mu_alpha}, Standard Deviation: {sigma_alpha}")

print(f"\nGaussian fit for C_reaction_Be_mass:")
print(f"Mean: {mu_Be}, Standard Deviation: {sigma_Be}")

# Plot the Gaussian fit on the existing reaction energy plots
x_vals = np.linspace(min(C_reaction_alpha_mass), max(C_reaction_alpha_mass), 1000)
y_vals_alpha = norm.pdf(x_vals, mu_alpha, sigma_alpha)

x_vals_Be = np.linspace(min(C_reaction_Be_mass), max(C_reaction_Be_mass), 1000)
y_vals_Be = norm.pdf(x_vals_Be, mu_Be, sigma_Be)

# Plot the results
plt.plot(C_reaction_alpha_mass, multiplied_hist_adjusted, label="12C reaction alpha energy")
plt.plot(x_vals, y_vals_alpha * np.max(multiplied_hist_adjusted), label="Gaussian Fit (alpha)", linestyle="--")

plt.plot(C_reaction_Be_mass, multiplied_hist_adjusted, label="12C reaction 8Be energy")
plt.plot(x_vals_Be, y_vals_Be * np.max(multiplied_hist_adjusted), label="Gaussian Fit (8Be)", linestyle="--")

plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('C12 reaction nuclei energies with Gaussian Fit')
plt.xlim(0, 2.5)
plt.show()





plt.plot(C_reaction_Be_mass, multiplied_hist_adjusted, label="12C reaction 8Be energy")

plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('C12 reaction nuclei energies with Gaussian Fit')
plt.xlim(0.7, 0.8)
plt.savefig('16O_energies.png', dpi=500)
plt.show()

plt.plot(C_reaction_alpha_mass, multiplied_hist_adjusted, label="12C reaction 8Be energy")

plt.legend()
plt.xlabel('Energy (MeV)')
plt.ylabel('Normalized Probability Density (multiplied by Beam Energy)')
plt.title('C12 reaction nuclei energies with Gaussian Fit')
plt.xlim(1.5, 1.6)
plt.savefig('C12_energies.png', dpi=500)
plt.show()
