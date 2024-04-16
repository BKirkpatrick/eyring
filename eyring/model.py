from scipy.constants import Planck, gas_constant, Boltzmann
from math import exp
import numpy as np
import json

gibbs = []  # Gibbs free energy
kappa = []  # transmission coefficient, [0 -> 1], determines fraction of flux through transition state which leads to product without recrossing transition state.
temperature = []  # Kelvin

def eyring_rate(temperature, gibbs, kappa):
    """Dissocation rate according to Eyring model"""
    k = (kappa * Boltzmann * temperature / Planck) * exp(- gibbs / gas_constant / temperature)
    return k

def calculate_gibbs(enthalpy, entropy, temperature):
    """Calculate the Gibbs free energy of system"""
    gibbs = enthalpy - (temperature * entropy)
    return gibbs

def extract_thermodynamic_parameters(seq, params_dict):
    """Calculate enthalpy and entropy of sequence"""
    # Convert the DNA sequence to its complementary sequence
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = seq.replace(" ", "")
    # Initialize total DeltaH and DeltaS
    dH = 0
    dS = 0
    # Iterate through the DNA sequence to get each nearest-neighbor pair
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]  # Original sequence pair

        if pair == "AA" or pair == "TT":
            key = "AA/TT"
        elif pair == "AT":
            key = "AT/AT"
        elif pair == "TA":
            key = "TA/TA"
        elif pair == "CA" or pair == "TG":
            key = "CA/TG"
        elif pair == "GT" or pair == "AC":
            key = "GT/AC"
        elif pair == "CT" or pair == "AG":
            key = "CT/AG"
        elif pair == "GA" or pair == "TC":
            key = "GA/TC"
        elif pair == "CG":
            key = "CG/CG"
        elif pair == "GC":
            key = "GC/GC"
        elif pair == "GG" or pair == "CC":
            key = "GG/CC"

    # add contribution from this pair
        dH += params_dict[key]['dH']
        dS += params_dict[key]['dS']
        
    return dH, dS

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
