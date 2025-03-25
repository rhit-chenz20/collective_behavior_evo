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
    
    neighbor_M = np.zeros((N, neighbor_size), dtype=int) # Initialize neighbor matrix for fixed size
    for i in range(N):
        nbs=[j for j in range(N) if i!=j]
        neighbors = np.random.choice(nbs, size=neighbor_size, replace=False) # Randomly select neighbors for interaction
        neighbor_M[i, :] = neighbors # Store neighbors for each individual
    
    # Simulation over generations
    for gen in range(generations):
        # Draw environmental effects (mean 0, variance env_std^2)
        e = np.random.normal(0, env_std, size=N)
        
        for i in range(N):
            nbs=[j for j in range(N) if i!=j] # List of individuals
            neighbors = neighbor_M[i,:] 
            for j in neighbors:
                ind1_a = sum(a[i,:])
                ind2_a = sum(a[j,:])
                ind1_e = e[i]
                ind2_e = e[j]
                z1 = (ind1_a + ind1_e + l * (ind2_a + ind2_e)) / (1 - l**2)
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
        # print(f"Generation {gen:3d}: Mean a = {mean_a:.3f}, Mean z = {mean_z:.3f}")
    
    return mean_a_list, mean_z_list

def run_replicates(N, generations, loci, l, mutation_std, beta, env_std, seed, replicates, neighbor_size):
    """
    Run multiple replicates of the simulation for the same parameter set.
    
    Parameters:
      replicates : Number of replicates to run.
    
    Returns:
      all_mean_a : List of lists containing mean additive genetic values for all replicates.
      all_mean_z : List of lists containing mean phenotypes for all replicates.
    """
    all_mean_a = []
    all_mean_z = []
    for rep in range(replicates):
        print(f"Running replicate {rep + 1}/{replicates}...")
        mean_a, mean_z = run_simulation(N, generations, loci, l, mutation_std, beta=beta, env_std=env_std, seed=seed, neighbor_size=neighbor_size)
        all_mean_a.append(mean_a)
        all_mean_z.append(mean_z)
    return all_mean_a, all_mean_z

if __name__ == "__main__":
    # Simulation parameters
    N = 20           # Population size 
    generations = 10  # Number of generations
    mutation_std = 0.05  # Mutation standard deviation for a
    beta = 0.1         # Selection gradient on phenotype
    env_std = 1.0      # Standard deviation of environmental noise
    seed = None        # Random seed for reproducibility
    loci = 100         # Number of loci for each individual
    replicates = 10    # Number of replicates to run
    neighbor_size = 10 # Fixed neighbor size for all simulations

    # Values of l to test
    l_values = [-0.75, -0.5, 0, 0.5, 0.75]

    # Dictionary to store results for each l value
    results = {}

    # Run simulations for each l value
    for l in l_values:
        print(f"Running simulations for l = {l}...")
        all_mean_a, all_mean_z = run_replicates(N, generations, loci, l, mutation_std, beta, env_std, seed, replicates, neighbor_size)
        results[l] = (all_mean_a, all_mean_z)

    # Define a color map for different l values
    colors = {
        -0.75: "blue",
        -0.5: "green",
        0: "black",
        0.5: "orange",
        0.75: "red"
    }

    # Plot the evolution of the mean additive genetic value and mean phenotype grouped by l values
    plt.figure(figsize=(12, 8))

    for l, (all_mean_a, all_mean_z) in results.items():
        color = colors[l]  # Get the color for the current l value
        # Plot individual replicates
        for rep in range(replicates):
            plt.plot(all_mean_a[rep], color=color, alpha=0.1, linestyle="--", label="_nolegend_")
            plt.plot(all_mean_z[rep], color=color, alpha=0.1, label="_nolegend_")
        # Plot mean across replicates
        mean_a = np.mean(all_mean_a, axis=0)
        mean_z = np.mean(all_mean_z, axis=0)
        plt.plot(mean_a, color=color, linestyle="--", label=f"Mean a (Psi={l})", linewidth=2)
        plt.plot(mean_z, color=color, label=f"Mean z (Psi={l})", linewidth=2)

    # Finalize plot
    plt.xlabel("Generation")
    plt.ylabel("Mean value")
    plt.title("Evolution of Interacting Phenotypes (Grouped by Psi values)")
    plt.legend()
    plt.tight_layout()
    plt.show()
