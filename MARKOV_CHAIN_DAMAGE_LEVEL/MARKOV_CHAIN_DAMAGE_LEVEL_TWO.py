###########################################################################################################
#                                         IN THE NAME OF ALLAH                                            #
#                STRUCTURAL DAMAGE RISK ASSESSMENT USING OPENSEES AND MARKOV CHAIN MODELING               #
#---------------------------------------------------------------------------------------------------------#
# ASSESSING STRUCTURAL DAMAGE RISKS USING A MARKOV CHAIN MODEL WITH THE OPENSEES LIBRARY                  #
# INVOLVES SIMULATING THE STRUCTURAL BEHAVIOR UNDER VARIOUS CONDITIONS AND THEN USING THE                 #
# RESULTS TO ESTIMATE TRANSITION PROBABILITIES BETWEEN DAMAGE STATES.                                     #
# BELOW IS AN EXAMPLE PYTHON CODE THAT DEMONSTRATES HOW TO USE OPENSEES                                   #
# TO SIMULATE A SIMPLE STRUCTURE AND THEN APPLY A MARKOV CHAIN MODEL TO ASSESS DAMAGE RISKS.              #
#---------------------------------------------------------------------------------------------------------#
# STEPS IN THE CODE:                                                                                      #
# [1] DEFINE THE STRUCTURE: CREATE A SIMPLE STRUCTURAL MODEL USING OPENSEES.                              #
# [2] SIMULATE DAMAGE STATES: PERFORM ANALYSES TO SIMULATE THE STRUCTURE'S RESPONSE                       #
# UNDER DIFFERENT LOADING CONDITIONS.                                                                     #
# [3] ESTIMATE TRANSITION PROBABILITIES: USE THE SIMULATION RESULTS TO ESTIMATE                           #
# TRANSITION PROBABILITIES BETWEEN DAMAGE STATES.                                                         #
# [4] APPLY MARKOV CHAIN MODEL: USE THE TRANSITION PROBABILITIES TO PREDICT FUTURE DAMAGE STATES.         #
#---------------------------------------------------------------------------------------------------------#
# EXPLANATION OF THE CODE:                                                                                #
# OPENSEES SIMULATION:                                                                                    #
# A SIMPLE 2D BEAM-COLUMN STRUCTURE IS MODELED.                                                           #
# THE STRUCTURE IS SUBJECTED TO A RANDOM LOAD LEVEL, AND THE DISPLACEMENT AT THE TOP NODE IS CALCULATED.  #
#---------------------------------------------------------------------------------------------------------#
# DAMAGE STATE CLASSIFICATION:                                                                            #
# [1] THE DISPLACEMENT IS USED TO CLASSIFY THE STRUCTURE INTO ONE OF THE DAMAGE STATES.                   #
# [2] TRANSITION PROBABILITY ESTIMATION:                                                                  #
# [3] THE STRUCTURE IS SIMULATED MULTIPLE TIMES UNDER RANDOM LOAD LEVELS.                                 #
# [4] THE TRANSITIONS BETWEEN DAMAGE STATES ARE COUNTED AND NORMALIZED TO                                 #
# ESTIMATE THE TRANSITION PROBABILITY MATRIX.                                                             #
#                                                                                                         #
# MARKOV CHAIN PREDICTION:                                                                                #
# THE ESTIMATED TRANSITION PROBABILITY MATRIX IS USED TO PREDICT THE FUTURE DAMAGE STATE PROBABILITIES.   #
#---------------------------------------------------------------------------------------------------------#
#                         THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                   #
#                                   EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                                   #
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
def simulate_structure(load_level):
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

    # Get displacement at node 2
    displacement = ops.nodeDisp(2, 1)
    base_reaction = -ops.eleResponse(1, 'force')[0]
    
    return displacement, base_reaction
#---------------------------------------------------------------------------------
# Function to classify damage state based on displacement
def classify_damage_state(displacement):
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
        DIS, BAS = simulate_structure(load_level)
        #print(DIS, BAS)
        current_state = classify_damage_state(DIS)

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
# Estimate transition probabilities
num_simulations = 1000000
transition_probabilities = ESTIMATE_TRANSITION_PROBABILITIES(num_simulations)
print("Estimated Transition Probability Matrix:")
print(transition_probabilities)

# Plot the Markov chain
PLOT_MARKOV_CHAIN(transition_probabilities)

# Predict future damage states
initial_state = 0  # Start with "No Damage"
num_steps = 10000  # Predict for 10 time steps
future_states = PREDICT_FUTURE_STATES(initial_state, transition_probabilities, num_steps)
print("\nPredicted Damage State Probabilities after", num_steps, "steps:")
for i, prob in enumerate(future_states):
	print(f"{DAMAGE_STATES[i]}: {prob:.4f}")
#---------------------------------------------------------------------------------    
