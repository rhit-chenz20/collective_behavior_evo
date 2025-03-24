import numpy as np
import matplotlib.pyplot as plt
# This script simulates the evolution of mean phenotype under reciprocal interactions
# between two traits, using the extended breeder's equation.

# Define parameters for simulation
num_generations = 50
population_size = 100
heritability = 0.5  # Direct genetic effect proportion
mutation_variance = 0.01  # Small genetic variance to maintain diversity

# Define a range of interaction effect coefficients (ψ) to analyze
psi_values = np.linspace(-0.99, 0.99, 31)  # Avoiding ψ = ±1 to prevent divergence

# Selection gradient (assuming a moderate selection pressure)
selection_gradient = 0.2

# Function to compute phenotype evolution using the extended breeder's equation
def compute_phenotypic_evolution(psi):
    mean_phenotypes = []
    mean_phenotype = 0  # Start from a mean of 0
    genetic_variance = 1  # Assume initial genetic variance is 1

    for generation in range(num_generations):
        denominator = 1 - psi**2  # The feedback effect in the reciprocal interaction model
        if np.abs(denominator) < 1e-6:  # Avoid division by very small numbers
            break

        # Extended breeder's equation: Δz̄ = (1 / (1 - ψ²)) * Gβ
        delta_mean_phenotype = (1 / denominator) * (genetic_variance * selection_gradient)

        # Update mean phenotype
        mean_phenotype += delta_mean_phenotype
        mean_phenotypes.append(mean_phenotype)

    return mean_phenotypes

# Store results for different ψ values
phenotypic_trajectories = {psi: compute_phenotypic_evolution(psi) for psi in psi_values}

# Plot the evolution of mean phenotype for different ψ values
plt.figure(figsize=(10, 6))

for psi, trajectory in phenotypic_trajectories.items():
    plt.plot(range(len(trajectory)), trajectory, label=f"ψ = {psi:.2f}")

plt.xlabel("Generation")
plt.ylabel("Mean Phenotypic Value")
plt.title("Evolution of Mean Phenotype under Reciprocal Interactions")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()
