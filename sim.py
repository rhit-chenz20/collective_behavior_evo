import numpy as np
import matplotlib.pyplot as plt


def run_simulation(N=20, generations=100, loci=10, l=0.75, mutation_std=0.05, beta=0.1, env_std=1.0, seed=None, neighbor_size=5):
    """
    Run the agent‐based simulation.
    
    Parameters:
      N             : Number of individuals (should be even).
      generations   : Number of generations to simulate.
      l             : Interaction effect coefficient (|l| < 1).
      mutation_std  : Standard deviation of mutation noise in a.
      beta          : Selection gradient (strength of selection on phenotype z).
      env_std       : Standard deviation of environmental effect.
      seed          : Random seed for reproducibility.
      
    Returns:
      mean_a_list   : List of mean additive genetic values (a) over generations.
      mean_z_list   : List of mean phenotypes (z) over generations.
    """
    if seed is not None:
        np.random.seed(seed)
        
    # Initialize the additive genetic values of 10 loci for N individuals (mean 0, variance = 1)
    a = np.random.normal(loc=0, scale=1, size=(N, loci))
    
    # initialize the phenotype matrix for N individuals (mean 0, variance = 1)
    z = np.zeros(((generations, N)))
    
    # Lists to record mean values over generations
    mean_a_list = []
    mean_z_list = []
    
    # Simulation over generations
    for gen in range(1,generations):
        # Draw environmental effects (mean 0, variance env_std^2)
        e = np.random.normal(0, env_std, size=N)
        
        for i in range(N):
            nbs=[j for j in range(N) if i!=j] # List of individuals
            neighbors = np.random.choice(nbs, size=neighbor_size) # Randomly select neighbors for interaction
            for j in neighbors:
                ind1_a = sum(a[i,:])
                ind2_a = sum(a[j,:])
                ind1_e = e[i]
                ind2_e = e[j]
                z1p = z[gen-1, j]
                z1 = (ind1_a + ind1_e + l * (ind2_a + ind2_e)) / (1 - l**2)
                # z1 = ind1_a + ind1_e + l * z1p
                z[gen, i] += z1
            z[gen, i] /= len(neighbors)  # Average over all other individuals
        
        # Compute fitness using an exponential function (ensuring positive fitness)
        fitness = np.exp(beta * z[gen, :])
        
        # Reproduce: sample N individuals with probability proportional to fitness.
        # Here we assume clonal reproduction: offspring inherit parent's a with mutation.
        parent_indices = np.random.choice(N, size=N, replace=True, p=fitness/fitness.sum())
        a_offspring = np.zeros((N, loci))
        # Note: mutation is applied to each locus independently
        for i in parent_indices:
            a_offspring[i, :] = a[i, :]+ np.random.normal(0, mutation_std, size=loci)
        # Update population genetic values
        a = a_offspring
        
        # Record means (note: mean environmental effect is zero so mean(z) ≈ mean(a)/(1-l))
        mean_a = np.mean([sum(asum) for asum in a])
        mean_z = np.mean(z[gen, :])
        mean_a_list.append(mean_a)
        mean_z_list.append(mean_z)
        
        # Optionally print progress
        print(f"Generation {gen:3d}: Mean a = {mean_a:.3f}, Mean z = {mean_z:.3f}")
    
    return mean_a_list, mean_z_list

if __name__ == "__main__":
    # Simulation parameters
    N = 20           # Population size (must be even)
    generations = 100  # Number of generations
    l = 0.75           # Interaction effect coefficient (e.g., 0.75 amplifies the response)
    mutation_std = 0.05  # Mutation standard deviation for a
    beta = 0.1         # Selection gradient on phenotype
    env_std = 1.0      # Standard deviation of environmental noise
    # seed = 42          # Random seed for reproducibility
    seed = None
    loci = 20          # Number of loci for each individual
    
    # Run the simulation
    mean_a, mean_z = run_simulation(N, generations, loci, l, mutation_std, beta, env_std, seed)

    # Plot the evolution of the mean additive genetic value and mean phenotype
    plt.figure(figsize=(8, 5))
    plt.plot(mean_a, label="Mean genetic value (a)")
    plt.plot(mean_z, label="Mean phenotype (z)")
    plt.xlabel("Generation")
    plt.ylabel("Mean value")
    plt.title("Evolution of Interacting Phenotypes")
    plt.legend()
    plt.tight_layout()
    plt.show()
