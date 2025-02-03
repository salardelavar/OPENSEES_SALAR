"""
Markov Chain Damage Level Simulation using OpenSeesPy and Bayesian Updating

This module simulates a simple 2D beam-column (representing a column)
in OpenSees, classifies its damage state based on displacement, and updates
a transition probability matrix via a Markov chain model. It then predicts
future damage states and plots the Markov chain.

Author: Salar Delavar Ghashghaei (Qashqai)
Email: salar.d.ghashghaei@gmail.com
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import logging

# Set up basic logging configuration
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


class MarkovChainDamageModel:
    DAMAGE_STATES = ["No Damage", "Minor Damage", "Moderate Damage", "Severe Damage", "Collapse"]

    def __init__(self, transition_matrix=None):
        # Default transition matrix (initial guess)
        if transition_matrix is None:
            self.transition_matrix = np.array([
                [0.55, 0.10, 0.20, 0.05, 0.10],  # No Damage
                [0.20, 0.65, 0.12, 0.00, 0.03],  # Minor Damage
                [0.15, 0.10, 0.50, 0.20, 0.05],  # Moderate Damage
                [0.10, 0.05, 0.10, 0.55, 0.20],  # Severe Damage
                [0.41, 0.50, 0.00, 0.00, 0.09]   # Collapse
            ])
        else:
            self.transition_matrix = transition_matrix

        # Structural and material properties (SI units)
        self.B = 0.40          # Width [m]
        self.H = 0.50          # Height [m]
        self.LENGTH = 6.0      # Column length [m]
        self.fc = 35           # Concrete compressive strength [N/mm²]
        self.Ec = 4700 * (self.fc)**0.5 * 1e6   # Young's modulus [Pa]
        self.Ar = self.B * self.H              # Cross-sectional area [m²]
        self.Ie = self.B * self.H**3 / 12       # Moment of inertia [m⁴]
        self.DENSITY = 2500    # Density [kg/m³]
        self.ele_mass = self.DENSITY * self.Ar  # Mass per unit length [kg/m]

    def simulate_structure(self, load_level):
        """
        Builds and analyzes a simple 2D beam-column model under a given load level.
        
        Parameters:
            load_level (float): The horizontal load applied at the top node [N].
        
        Returns:
            tuple: (displacement, base_reaction) where displacement is the horizontal
                   displacement at node 2 [m] and base_reaction is the reaction force at the base [N].
        """
        # Clear previous model
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # Define nodes: node1 (base) fixed, node2 (top) free
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, self.LENGTH)

        # Set boundary conditions: fix all DOFs at node 1
        ops.fix(1, 1, 1, 1)

        # Define element: elasticBeamColumn element with geometric nonlinearity
        ops.geomTransf('Linear', 1)
        # Create an elastic beam-column element with id 1
        ops.element('elasticBeamColumn', 1, 1, 2, self.Ar, self.Ec, self.Ie, 1, '-mass', self.ele_mass)

        # Apply load pattern: horizontal load at node 2
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(2, load_level, 0, 0)

        # Analysis settings
        ops.system('BandGeneral')
        ops.numberer('Plain')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1)
        ops.algorithm('Newton')
        ops.analysis('Static')
        ops.analyze(1)

        displacement = ops.nodeDisp(2, 1)  # horizontal displacement at node 2
        # Base reaction extracted from the element force (assuming forceBeamColumn returns [Px, Py, M])
        base_reaction = -ops.eleResponse(1, 'force')[0]

        return displacement, base_reaction

    @staticmethod
    def classify_damage_state(displacement):
        """
        Classify the damage state based on displacement magnitude.
        
        Parameters:
            displacement (float): The horizontal displacement [m].
        
        Returns:
            int: Index representing the damage state.
        """
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

    def bayesian_update_transition_probabilities(self, num_simulations, prior_counts=None):
        """
        Runs multiple simulations to update the transition probability matrix using Bayesian updating.
        
        Parameters:
            num_simulations (int): Number of simulation runs.
            prior_counts (np.ndarray): Prior counts for the transition matrix. If None, uniform prior is used.
        
        Returns:
            np.ndarray: The updated transition probability matrix.
        """
        num_states = len(self.DAMAGE_STATES)
        if prior_counts is None:
            transition_counts = np.ones((num_states, num_states))  # Uniform Dirichlet prior
        else:
            transition_counts = prior_counts.copy().astype(float)

        for _ in range(num_simulations):
            load_level = np.random.uniform(0, 3)  # Random load level between 0 and 3 N (or scaled units)
            displacement, _ = self.simulate_structure(load_level)
            current_state = self.classify_damage_state(displacement)
            # Simulate next state based on the assumed true transition probabilities
            next_state = np.random.choice(num_states, p=self.transition_matrix[current_state])
            transition_counts[current_state, next_state] += 1

        # Normalize counts to obtain probabilities
        updated_probabilities = transition_counts / transition_counts.sum(axis=1, keepdims=True)
        return updated_probabilities

    @staticmethod
    def predict_future_states(initial_state, transition_matrix, num_steps):
        """
        Predicts the future damage state probabilities after a number of steps using the Markov chain.
        
        Parameters:
            initial_state (int): The index of the initial damage state.
            transition_matrix (np.ndarray): The transition probability matrix.
            num_steps (int): Number of Markov chain steps.
        
        Returns:
            np.ndarray: A probability distribution over the damage states.
        """
        state_prob = np.zeros(len(MarkovChainDamageModel.DAMAGE_STATES))
        state_prob[initial_state] = 1.0
        for _ in range(num_steps):
            state_prob = np.dot(state_prob, transition_matrix)
        return state_prob

    @staticmethod
    def plot_markov_chain(transition_matrix, damage_states):
        """
        Plots the Markov chain as a directed graph with edge labels showing transition probabilities.
        
        Parameters:
            transition_matrix (np.ndarray): Transition probability matrix.
            damage_states (list): List of state names.
        """
        G = nx.DiGraph()
        # Add nodes
        for state in damage_states:
            G.add_node(state)
        # Add edges with nonzero probabilities
        num_states = len(damage_states)
        for i in range(num_states):
            for j in range(num_states):
                if transition_matrix[i, j] > 0:
                    G.add_edge(damage_states[i], damage_states[j], weight=transition_matrix[i, j])
        pos = nx.spring_layout(G)
        plt.figure(figsize=(10, 6))
        nx.draw(G, pos, with_labels=True, node_size=3000, node_color='lightblue', font_size=10, font_weight='bold')
        edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
        plt.title("Markov Chain Model of Structural Damage States")
        plt.show()


def main():
    model = MarkovChainDamageModel()
    num_simulations = 100000  # You can adjust the number of simulations as needed

    logging.info("Starting Bayesian update of transition probabilities...")
    updated_transitions = model.bayesian_update_transition_probabilities(num_simulations)
    logging.info("Updated Transition Probability Matrix:")
    print(updated_transitions)

    # Plot the updated Markov chain
    model.plot_markov_chain(updated_transitions, model.DAMAGE_STATES)

    # Predict future damage state probabilities
    initial_state = 0  # Start with "No Damage"
    num_steps = 10000
    future_state_prob = model.predict_future_states(initial_state, updated_transitions, num_steps)
    logging.info("Predicted Damage State Probabilities after Bayesian Update:")
    for state, prob in zip(model.DAMAGE_STATES, future_state_prob):
        print(f"{state}: {prob:.4f}")


if __name__ == '__main__':
    main()
