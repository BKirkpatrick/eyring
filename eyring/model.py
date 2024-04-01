from scipy.constants import Planck, gas_constant, Boltzmann
from math import exp
import numpy as np

gibbs = []  # Gibbs free energy
kappa = []  # transmission coefficient, [0 -> 1], determines fraction of flux through transition state which leads to product without recrossing transition state.
temperature = []  # Kelvin

def eyring_rate(temperature, gibbs, kappa):
    """Dissocation rate according to Eyring model"""
    k = (kappa * Boltzmann * temperature / Planck) * exp(- gibbs / gas_constant / temperature)
    return k

def calculate_gibbs(enthalpy, entropy, temperature):
    """Calculate the Gibbs free energy of system"""
    gibbs = enthalpy - temperature * entropy
    return gibbs

def extract_thermodynamic_parameters(dna_sequence):
    """Calculate enthalpy and entropy of sequence"""
    # Thermodynamic parameters for nearest-neighbor pairs
    parameters = {
        "AA/TT": {"DeltaH": 7.9, "DeltaS": 22.2},
        "AT/TA": {"DeltaH": 7.2, "DeltaS": 20.4},
        "TA/AT": {"DeltaH": 7.2, "DeltaS": 21.3},
        "CA/GT": {"DeltaH": 8.5, "DeltaS": 22.7},
        "GT/CA": {"DeltaH": 8.4, "DeltaS": 22.4},
        "CT/GA": {"DeltaH": 7.8, "DeltaS": 21.0},
        "GA/CT": {"DeltaH": 8.2, "DeltaS": 22.2},
        "CG/GC": {"DeltaH": 10.6, "DeltaS": 27.2},
        "GC/CG": {"DeltaH": 9.8, "DeltaS": 24.4},
        "GG/CC": {"DeltaH": 8.0, "DeltaS": 19.9}
    }

    # Convert the DNA sequence to its complementary sequence
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    comp_sequence = ''.join([complement[base] for base in dna_sequence])

    # Initialize total DeltaH and DeltaS
    total_DeltaH = 0
    total_DeltaS = 0

    # Iterate through the DNA sequence to get each nearest-neighbor pair
    for i in range(len(dna_sequence) - 1):
        pair = dna_sequence[i:i+2]  # Original sequence pair
        comp_pair = comp_sequence[i:i+2]  # Complementary sequence pair
        nn_pair = pair + '/' + comp_pair  # Construct the nearest-neighbor key

        # Look up the thermodynamic parameters for the current nearest-neighbor pair
        if nn_pair in parameters:
            total_DeltaH += 1000 * parameters[nn_pair]["DeltaH"]
            total_DeltaS += parameters[nn_pair]["DeltaS"]
        elif nn_pair[::-1] in parameters:
            total_DeltaH += 1000 * parameters[nn_pair[::-1]]["DeltaH"]
            total_DeltaS += parameters[nn_pair[::-1]]["DeltaS"]
        else:
            print(f"Parameters for {nn_pair} not found.")

    return total_DeltaH, total_DeltaS

def simulate_dna_melting(num_molecules, rate_constant, num_steps):
    # Initialize all molecules as not melted (0)
    states = np.ones(num_molecules)
    
    # Simulation loop
    for step in range(num_steps):
        for i in range(num_molecules):
            if states[i] == 1:  # If not already melted
                if np.random.rand() < rate_constant:  # Compare to rate constant
                    states[i] = 0  # Melt the DNA molecule
                    
    return states
