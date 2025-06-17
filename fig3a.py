import numpy as np
import matplotlib.pyplot as plt

# reproduce figure 1(2D)

# Define parameters
genetic_variance = 2  # Assume initial genetic variance is 1
selection_gradient = 1.0  # Moderate selection pressure
psi_values = np.linspace(-0.9, 0.9, 50)  # Range of ψ values excluding ±1

# Compute expected change in phenotype in one generation using the extended breeder's equation
delta_mean_phenotype = (1 / (1 - psi_values**2))**2

# Plot the expected change in one generation as a function of ψ
plt.figure(figsize=(8, 5))
plt.plot(psi_values, delta_mean_phenotype, marker='o', linestyle='-')
plt.xlabel("Interaction Effect Coefficient (ψ11)")
plt.ylabel("Relative Change in Mean Phenotype")
plt.title("Expected Phenotypic Change in One Generation")
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Reference line at Δz̄ = 0
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Reference line at ψ = 0
plt.grid(True)
plt.show()