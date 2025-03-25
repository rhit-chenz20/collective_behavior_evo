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
    a1 = np.random.normal(loc=0, scale=1, size=(N, loci))
    a2 = np.random.normal(loc=0, scale=1, size=(N, loci))
    
    # initialize the phenotype matrix for N individuals (mean 0, variance = 1)
    z1 = np.zeros(((generations, N)))
    z2 = np.zeros(((generations, N)))
    
    
    # Lists to record mean values over generations
    mean_a1_list = []
    mean_a2_list = []
    mean_z1_list = []
    mean_z2_list = []
    
    # Simulation over generations
    for gen in range(generations):
        # Draw environmental effects (mean 0, variance env_std^2)
        e1 = np.random.normal(0, env_std, size=N)
        e2 = np.random.normal(0, env_std, size=N)
        
        for i in range(N):
            nbs=[j for j in range(N) if i!=j] # List of individuals
            neighbors = np.random.choice(nbs, size=neighbor_size) # Randomly select neighbors for interaction
            for j in neighbors:
                ind1_a1 = sum(a1[i,:])
                ind2_a1 = sum(a1[j,:])
                ind1_a2 = sum(a2[i,:])
                ind2_a2 = sum(a2[j,:])
                ind1_e1 = e1[i]  
                ind2_e1 = e1[j]
                ind1_e2 = e2[i]
                ind2_e2 = e2[j]
                
                z1f = ind1_a1 + ind1_e1 + l * (ind2_a2 + ind2_e2)
                z1[gen, i] += z1f
            z2f = ind1_a2 + ind1_e2
            z2[gen, i] = z2f
            z1[gen, i] /= len(neighbors)  # Average over all other individuals
            
        # Compute fitness using an exponential function (ensuring positive fitness)
        fitness1 = np.exp(beta * z1[gen, :])
        fitness2 = np.exp(beta * z2[gen, :])
        
        fitness = (fitness1 + fitness2) / 2
        
        # Reproduce: sample N individuals with probability proportional to fitness.
        # Here we assume clonal reproduction: offspring inherit parent's a with mutation.
        parent_indices = np.random.choice(N, size=N, replace=True, p=fitness/fitness.sum())
        a1_offspring = np.zeros((N, loci))
        a2_offspring = np.zeros((N, loci))
        # Note: mutation is applied to each locus independently
        for i in parent_indices:
            a1_offspring[i, :] = a1[i, :]+ np.random.normal(0, mutation_std, size=loci)
            a2_offspring[i, :] = a2[i, :]+ np.random.normal(0, mutation_std, size=loci)
        # Update population genetic values
        a1 = a1_offspring
        a2 = a2_offspring
        
        # Record means (note: mean environmental effect is zero so mean(z) ≈ mean(a)/(1-l))
        mean_a1 = np.mean([sum(asum) for asum in a1])
        mean_a2 = np.mean([sum(asum) for asum in a2])
        mean_z1 = np.mean(z1[gen, :])
        mean_z2 = np.mean(z2[gen, :])
        mean_a1_list.append(mean_a1)
        mean_a2_list.append(mean_a2)
        mean_z1_list.append(mean_z1)
        mean_z2_list.append(mean_z2)
        
        # Optionally print progress
        # print(f"Generation {gen:3d}: Mean a = {mean_a:.3f}, Mean z = {mean_z:.3f}")
    
    return mean_a1_list,mean_a2_list, mean_z1_list, mean_z2_list

def run_replicates(N, generations, loci, l, mutation_std, beta, env_std, seed, replicates, neighbor_size):
    """
    Run multiple replicates of the simulation for the same parameter set.
    
    Parameters:
      replicates : Number of replicates to run.
    
    Returns:
      all_mean_a1 : List of lists containing mean additive genetic values for trait 1 across all replicates.
      all_mean_a2 : List of lists containing mean additive genetic values for trait 2 across all replicates.
      all_mean_z1 : List of lists containing mean phenotypes for trait 1 across all replicates.
      all_mean_z2 : List of lists containing mean phenotypes for trait 2 across all replicates.
    """
    all_mean_a1 = []
    all_mean_a2 = []
    all_mean_z1 = []
    all_mean_z2 = []
    for rep in range(replicates):
        print(f"Running replicate {rep + 1}/{replicates}...")
        mean_a1, mean_a2, mean_z1, mean_z2 = run_simulation(N, generations, loci, l, mutation_std, beta, env_std, seed, neighbor_size)
        all_mean_a1.append(mean_a1)
        all_mean_a2.append(mean_a2)
        all_mean_z1.append(mean_z1)
        all_mean_z2.append(mean_z2)
    return all_mean_a1, all_mean_a2, all_mean_z1, all_mean_z2

if __name__ == "__main__":
    # Simulation parameters
    N = 20           # Population size 
    generations = 30  # Number of generations
    mutation_std = 0.05  # Mutation standard deviation for a
    beta = 0.1         # Selection gradient on phenotype
    env_std = 1.0      # Standard deviation of environmental noise
    seed = None        # Random seed for reproducibility
    loci = 100         # Number of loci for each individual
    replicates = 30    # Number of replicates to run
    neighbor_size = 10 # Fixed neighbor size for all simulations

    # Values of l to test
    l_values = [-0.75, -0.5, 0, 0.5, 0.75]

    # Dictionary to store results for each l value
    results = {}

    # Run simulations for each l value
    for l in l_values:
        print(f"Running simulations for l = {l}...")
        all_mean_a1, all_mean_a2, all_mean_z1, all_mean_z2 = run_replicates(
            N, generations, loci, l, mutation_std, beta, env_std, seed, replicates, neighbor_size
        )
        results[l] = (all_mean_a1, all_mean_a2, all_mean_z1, all_mean_z2)

    # Define a color map for different l values
    colors = {
        -0.75: "blue",
        -0.5: "green",
        0: "black",
        0.5: "orange",
        0.75: "red"
    }

    # Plot the evolution of the mean additive genetic values and mean phenotypes grouped by l values
    plt.figure(figsize=(12, 8))

    for l, (all_mean_a1, all_mean_a2, all_mean_z1, all_mean_z2) in results.items():
        color = colors[l]  # Get the color for the current l value
        # Plot individual replicates for trait 1
        for rep in range(replicates):
            plt.plot(all_mean_a1[rep], color=color, alpha=0.1, linestyle="--", label="_nolegend_")
            plt.plot(all_mean_z1[rep], color=color, alpha=0.1, label="_nolegend_")
        # Plot mean across replicates for trait 1
        mean_a1 = np.mean(all_mean_a1, axis=0)
        mean_z1 = np.mean(all_mean_z1, axis=0)
        plt.plot(mean_a1, color=color, linestyle="--", label=f"Mean a1 (Psi={l})", linewidth=2)
        plt.plot(mean_z1, color=color, label=f"Mean z1 (Psi={l})", linewidth=2)

    # Finalize plot for trait 1
    plt.xlabel("Generation")
    plt.ylabel("Mean value")
    plt.title("Evolution of Interacting Phenotypes (Trait 1, Grouped by Psi values)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot the evolution of the mean additive genetic values and mean phenotypes for trait 2
    plt.figure(figsize=(12, 8))

    for l, (all_mean_a1, all_mean_a2, all_mean_z1, all_mean_z2) in results.items():
        color = colors[l]  # Get the color for the current l value
        # Plot individual replicates for trait 2
        for rep in range(replicates):
            plt.plot(all_mean_a2[rep], color=color, alpha=0.1, linestyle="--", label="_nolegend_")
            plt.plot(all_mean_z2[rep], color=color, alpha=0.1, label="_nolegend_")
        # Plot mean across replicates for trait 2
        mean_a2 = np.mean(all_mean_a2, axis=0)
        mean_z2 = np.mean(all_mean_z2, axis=0)
        plt.plot(mean_a2, color=color, linestyle="--", label=f"Mean a2 (Psi={l})", linewidth=2)
        plt.plot(mean_z2, color=color, label=f"Mean z2 (Psi={l})", linewidth=2)

    # Finalize plot for trait 2
    plt.xlabel("Generation")
    plt.ylabel("Mean value")
    plt.title("Evolution of Interacting Phenotypes (Trait 2, Grouped by Psi values)")
    plt.legend()
    plt.tight_layout()
    plt.show()
