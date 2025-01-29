###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#---------------------------------------------------------------------------------------------------------#
#                    STRUCTURAL DAMAGE RISK ASSESSMENT USING OPENSEES AND MARKOV CHAIN MODELING           #
#---------------------------------------------------------------------------------------------------------#
# This program assesses structural damage risks using a Markov Chain Model integrated with OpenSees.      #
# The simulation evaluates structural behavior under varying conditions and estimates transition          #
# probabilities between different damage states. The key steps of the process are outlined below:         #
#---------------------------------------------------------------------------------------------------------#
#                                        PROGRAM WORKFLOW                                                 #
#---------------------------------------------------------------------------------------------------------#
# [1] STRUCTURAL MODEL DEFINITION:                                                                        #
#    - A simple 2D beam-column model is created in OpenSees.                                              #
#    - Material properties and geometric parameters are defined.                                          #
#    - Boundary conditions and loading conditions are applied.                                            #
#---------------------------------------------------------------------------------------------------------#
# [2] DAMAGE STATE SIMULATION:                                                                            #
#    - The structure undergoes multiple analyses under varying load levels.                               #
#    - The displacement at the top node is recorded.                                                      #
#    - Based on displacement values, the structure is classified into one of five damage states:          #
#      [1] No Damage   [2] Minor Damage   [3] Moderate Damage   [4] Severe Damage   [5] Collapse          #
#---------------------------------------------------------------------------------------------------------#
# [3] TRANSITION PROBABILITY ESTIMATION:                                                                  #
#    - Multiple simulations are run to track state transitions.                                           #
#    - Transition frequencies are recorded and normalized to form the transition probability matrix.      #
#---------------------------------------------------------------------------------------------------------#
# [4] MARKOV CHAIN PREDICTION:                                                                            #
#    - The estimated transition probability matrix is used to forecast the probability of future states.  #
#    - Bayesian updating is applied to refine transition probabilities based on new observations.         #
#---------------------------------------------------------------------------------------------------------#
#                                        IMPLEMENTATION DETAILS                                           #
#---------------------------------------------------------------------------------------------------------#
# - OpenSees is used for structural simulation.                                                           #
# - NumPy and NetworkX are employed for numerical computations and graph visualization.                   #
# - The Markov Chain Model predicts future structural conditions based on observed transitions.           #
# - The Bayesian approach updates transition probabilities dynamically as more data becomes available.    #
#---------------------------------------------------------------------------------------------------------#
#                         THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                   #
#                                EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                                      #
###########################################################################################################
#--------------------------------------------------------------------
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
#--------------------------------------------------------------------
# Define parameters (units: m, N)
# Define material properties
# Define a rectangular concrete section using fibers.
B = 0.40                     # [m] Width of the rectangular section
H = 0.50                     # [m] Height of the rectangular section
LENGTH = 6                   # [m] Column length
fc = 35                      # [N/mm²] Concrete Compressive strength
Ec = 4700 * (fc)**0.5 * 1e6  # [Pa] Young's Modulus for concrete
Ar = B * H                   # [m^2] Cross-sectional area
Ie = B * H**3 /12            # [m^4] Moment of Inertia
DENSITY = 2500               # [kg/m²] Concrete Density
Ele_Mass = DENSITY * Ar      # [kg/m] Mass per unit length
#---------------------------------------------------------------------------------
# Define the Markov chain states
DAMAGE_STATES = ["No Damage", "Minor Damage", "Moderate Damage", "Severe Damage", "Collapse"]
#---------------------------------------------------------------------------------
# Define the transition probability matrix (initial guess)
TRANSITION_MATRIX = np.array([
    [0.55, 0.10, 0.20, 0.05, 0.10],  # No Damage
    [0.20, 0.65, 0.12, 0.00, 0.03],  # Minor Damage
    [0.15, 0.10, 0.50, 0.20, 0.05],  # Moderate Damage
    [0.10, 0.05, 0.10, 0.55, 0.20],  # Severe Damage
    [0.41, 0.50, 0.00, 0.00, 0.09]   # Collapse
])
#---------------------------------------------------------------------------------
# Function to simulate the structure using OpenSees
def SIMULATE_STRUCTURE(load_level):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # Define nodes
    node1 = [0, 0]
    node2 = [0, LENGTH]

    ops.node(1, *node1)
    ops.node(2, *node2)

    # Define elements
    element_id = 1
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', element_id, 1, 2, Ar, Ec, Ie, 1,'-mass', Ele_Mass)

    # Define boundary conditions
    ops.fix(1, 1, 1, 1)

    # Apply load
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, load_level, 0, 0)

    # Perform analysis
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1)
    ops.algorithm('Newton')
    ops.analysis('Static')
    ops.analyze(1)

    displacement = ops.nodeDisp(2, 1)                 # Displacement at node 2
    base_reaction = -ops.eleResponse(1, 'force')[0]   # Base-reaction at node 1
    
    return displacement, base_reaction
#---------------------------------------------------------------------------------
# Function to classify damage state based on displacement
def CLASSIFY_DAMAGE_STATE(displacement):
    if displacement < 0.1:
        return 0  # No Damage
    elif displacement < 2.5:
        return 1  # Minor Damage
    elif displacement < 5.1:
        return 2  # Moderate Damage
    elif displacement < 8.3:
        return 3  # Severe Damage
    else:
        return 4  # Collapse
#---------------------------------------------------------------------------------
# Simulate the structure under different load levels
def ESTIMATE_TRANSITION_PROBABILITIES(num_simulations):
    transition_counts = np.zeros((len(DAMAGE_STATES), len(DAMAGE_STATES)))

    for _ in range(num_simulations):
        # Random load level between 0 and 3
        load_level = np.random.uniform(0, 3)
        DIS, BAS = SIMULATE_STRUCTURE(load_level)
        #print(DIS, BAS)
        current_state = CLASSIFY_DAMAGE_STATE(DIS)

        # Simulate next state based on current state
        next_state = np.random.choice(len(DAMAGE_STATES), p=TRANSITION_MATRIX[current_state])
        transition_counts[current_state, next_state] += 1

    # Normalize transition counts to probabilities
    transition_probabilities = transition_counts / transition_counts.sum(axis=1, keepdims=True)
    return transition_probabilities
#---------------------------------------------------------------------------------
# Function to predict future damage states using Markov chain
def PREDICT_FUTURE_STATES(initial_state, transition_matrix, num_steps):
    state_probabilities = np.zeros(len(DAMAGE_STATES))
    state_probabilities[initial_state] = 1.0

    for _ in range(num_steps):
        state_probabilities = np.dot(state_probabilities, transition_matrix)

    return state_probabilities
#---------------------------------------------------------------------------------
# Function to plot the Markov chain
def PLOT_MARKOV_CHAIN(transition_matrix):
    G = nx.DiGraph()

    # Add nodes
    for i, state in enumerate(DAMAGE_STATES):
        G.add_node(state)

    # Add edges with transition probabilities
    for i in range(len(DAMAGE_STATES)):
        for j in range(len(DAMAGE_STATES)):
            if transition_matrix[i, j] > 0:
                G.add_edge(DAMAGE_STATES[i], DAMAGE_STATES[j], weight=transition_matrix[i, j])

    # Plot the graph
    pos = nx.spring_layout(G)
    plt.figure(figsize=(10, 6))
    nx.draw(G, pos, with_labels=True, node_size=3000, node_color='lightblue', font_size=10, font_weight='bold')
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
    plt.title("Markov Chain Model of Structural Damage States")
    plt.show()
#---------------------------------------------------------------------------------
# Add this function for Bayesian updating
def BAYESIAN_UPDATE_TRANSITION_PROBABILITIES(num_simulations, prior_counts):
    # Initialize transition counts with prior parameters
    transition_counts = prior_counts.copy().astype(float)
    
    for _ in range(num_simulations):
        # Random load level between 0 and 3
        load_level = np.random.uniform(0, 3)
        DIS, BAS = SIMULATE_STRUCTURE(load_level)
        current_state = CLASSIFY_DAMAGE_STATE(DIS)
        
        # Simulate next state based on true transition matrix (from original code)
        next_state = np.random.choice(len(DAMAGE_STATES), p=TRANSITION_MATRIX[current_state])
        
        # Update counts for Bayesian posterior
        transition_counts[current_state, next_state] += 1
    
    # Calculate posterior transition probabilities
    transition_probabilities = transition_counts / transition_counts.sum(axis=1, keepdims=True)
    return transition_probabilities
#---------------------------------------------------------------------------------
# Define Dirichlet prior parameters (uniform prior)
PRIOR_COUNTS = np.ones((len(DAMAGE_STATES), len(DAMAGE_STATES)))

# Estimate transition probabilities with Bayesian updating
num_simulations = 1000000
updated_transition_probs = BAYESIAN_UPDATE_TRANSITION_PROBABILITIES(num_simulations, PRIOR_COUNTS)

print("Bayesian Updated Transition Probability Matrix:")
print(updated_transition_probs)

# Plot the updated Markov chain
PLOT_MARKOV_CHAIN(updated_transition_probs)

# Predict future states using updated probabilities
initial_state = 0  # Start with "No Damage"
num_steps = 10000
future_states = PREDICT_FUTURE_STATES(initial_state, updated_transition_probs, num_steps)

print("\nPredicted Damage State Probabilities after Bayesian Update:")
for i, prob in enumerate(future_states):
    print(f"{DAMAGE_STATES[i]}: {prob:.4f}")
#---------------------------------------------------------------------------------    
