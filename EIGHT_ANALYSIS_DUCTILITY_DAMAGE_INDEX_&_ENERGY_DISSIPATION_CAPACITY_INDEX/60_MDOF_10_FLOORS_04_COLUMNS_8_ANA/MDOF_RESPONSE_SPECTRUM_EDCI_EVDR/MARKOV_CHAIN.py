"""
This Python script performs a structural damage analysis by evaluating displacement data from an input file
 (or directly supplied data) and visualizing the results using a Markov Chain model. Below is a breakdown of each part of the code:

1. Input Data:
The script starts by reading displacement data from a text file or using directly provided data (`DATA`). The displacement data represents the movement of a structure, typically caused by external forces such as seismic activity (e.g., earthquakes). This data is assumed to be in a single column of numerical values.
- `np.loadtxt(file_path)` is used to read the data from a specified file (e.g., `Ground_Acceleration_1.txt`), or if no file is specified (`FILE_TF == False`), the script will use data directly supplied through the `DATA` parameter.
- This data represents measurements of structural displacements over time, which can vary depending on environmental or operational factors.

2. Quantile-Based Thresholds:
The script divides the displacement data into different levels of structural damage based on specific quantiles (percentiles). This quantile-based approach dynamically adapts to the distribution of the input data, making it flexible for various types of displacement datasets.
- Q1 (5th Percentile): Represents the "No Damage" state. This is the lowest 5% of the displacement values, indicating minimal movement and no noticeable damage to the structure.
- Q2 (10th Percentile): Represents the "Minor Damage" state. This range indicates slight deformations that may not affect structural integrity significantly.
- Q3 (50th Percentile): Represents the "Moderate Damage" state. This threshold captures noticeable deformations that could require repair.
- Q4 (80th Percentile): Represents the "Severe Damage" state. This range marks major deformations, suggesting serious structural problems.
- Q5 (95th Percentile): Represents the "Failure" state, where extreme displacement values indicate that the structure has reached or exceeded its capacity and is at risk of collapse.

These thresholds are determined using the `np.quantile()` function, which divides the data into quantile-based intervals. These quantile levels create a dynamic way to assess the severity of structural displacement based on the actual data.

3. Markov Transition Matrix:
Once the data is divided into states, the script calculates a Markov Transition Matrix, which represents the probabilities of transitions between the states over time.
- Transition Matrix: A matrix where each element \( P[i, j] \) represents the probability of transitioning from state `i` to state `j`. For example, it might show the probability of a structure staying in the "No Damage" state or progressing to "Minor Damage".
- Counting Transitions: The script counts how many times the structure transitions from one state to another by iterating over the displacement data and checking how the state changes from one time step to the next.
- Normalization: After counting the transitions, the matrix is normalized (i.e., transformed into probabilities) by dividing each row by its total sum. This ensures that the probabilities for each state sum to 1.

4. NetworkX Visualization:
To better understand the Markov Transition Matrix, the script visualizes it as a directed graph using NetworkX.
- Nodes: Represent the different states (e.g., "No Damage", "Minor Damage", etc.).
- Edges: Represent transitions between states. The thickness of the edges corresponds to the transition probability, and each edge is labeled with the exact transition probability in percentage form. 
- Spring Layout: The positions of the nodes are arranged using a spring layout to minimize edge overlaps and make the graph more readable. 
- Edge Labels: The edges are labeled with the transition probabilities in percentage format (e.g., "15.00000 %"), providing insight into the likelihood of each transition.

Key Observations in the Graph:
- Dominant Transitions: Transitions that occur frequently (like "No Damage → No Damage") will have thicker edges.
- Escalation Trends: If there is a strong probability of transitioning from a lower-damage state (like "Minor Damage") to a higher-damage state (like "Moderate Damage"), the graph will highlight these edges.
- Recovery Trends: If the data indicates any recovery (e.g., transitioning from "Severe Damage" to "Moderate Damage"), these edges will also be visualized.



5. Practical Implications:
The output provides valuable insights for engineers and disaster response teams:
- Engineering Insight: By understanding how likely different damage states are, engineers can predict how the structure will behave under continuous stress or seismic activity. If a structure has a high probability of moving from "Severe Damage" to "Failure", it suggests that immediate intervention is needed. 
- Disaster Response: By analyzing the transition matrix and visualizing the critical thresholds, teams can prepare for rapid escalation in damage. This can assist in resource allocation and focus efforts on the most at-risk structures.
- Structural Design: If designers understand the transition probabilities, they can take steps to enhance resilience in areas where damage escalation is likely (e.g., increasing reinforcement in sections of a building that are more likely to experience higher displacement).



6. Limitations:
While this approach provides valuable insights, there are some inherent limitations:
- Static Data: The model assumes that transition probabilities remain constant over time, which may not be realistic in dynamic, real-world scenarios where factors such as fatigue, repair, or additional stresses can change the system's behavior. 
- Dataset Dependency: The results depend heavily on the quality and size of the input data (`Ground_Acceleration_1.txt`). A small or biased dataset may lead to misleading conclusions.
- Simplified States: The model simplifies a complex phenomenon (structural behavior) into discrete states. This simplification may overlook nuanced behaviors that occur in between these states, making the model less precise in capturing the full complexity of structural dynamics.

In summary, this script models the transition of a structure between different damage states under stress using a Markov chain and provides a visual representation of the probabilities for each state transition. This can help in structural assessment, disaster management, and design improvements, although it does make certain assumptions that may not always apply in real-world cases.
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def MARKOV_CHAIN(FILE_TF, file_path, DATA):
    # Load the displacement data from the provided text file
    if FILE_TF == True:
        data = np.loadtxt(file_path)  # Ensure the file contains a single column of values (displacement, velocity, acceleration, base-reaction)
    if FILE_TF == False:
        data = DATA    #  (displacement, velocity, acceleration, base-reaction)
    # Define state boundaries based on quantiles
    Q1 = np.quantile(data, 0.05)  # No Damage Level
    Q2 = np.quantile(data, 0.10)  # Minor Damage Level
    Q3 = np.quantile(data, 0.50)  # Moderate Damage Level
    Q4 = np.quantile(data, 0.80)  # Severe Damage Level
    #Q5 = np.quantile(data, 0.95)  # Failure Level

    # Define critical thresholds and labels for states
    state_limits = [Q1, Q2, Q3, Q4]  # Critical thresholds for damage levels
    state_labels = ["No Damage", "Minor Damage", "Moderate Damage", "Severe Damage", "Failure"]
    num_states = len(state_limits) + 1  # Total number of states (thresholds + 1)

    # Function to determine the state of a sample based on displacement
    def get_state(value, limits):
        for i, limit in enumerate(limits):
            if value <= limit:
                return i  # Return the corresponding state index
        return len(limits)  # If the value exceeds all limits, assign to the last state

    # Map each displacement value in the data to a state
    states = [get_state(disp, state_limits) for disp in data]

    # Initialize the Markov transition matrix with zeros
    transition_matrix = np.zeros((num_states, num_states))

    # Count transitions between states
    for i in range(len(states) - 1):
        transition_matrix[states[i], states[i + 1]] += 1

    # Normalize the transition matrix to represent probabilities
    row_sums = transition_matrix.sum(axis=1, keepdims=True)  # Calculate the sum of each row
    transition_matrix = np.divide(transition_matrix, row_sums, where=row_sums != 0)  # Normalize rows where sum is not zero

    # Print the Markov Transition Matrix
    print("Markov Transition Matrix:")
    print(transition_matrix)

    # Visualize the transition matrix using NetworkX
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes and edges with weights based on transition probabilities
    for i in range(num_states):
        for j in range(num_states):
            if transition_matrix[i, j] > 0:  # Add edges only if the probability is non-zero
                G.add_edge(state_labels[i], state_labels[j], weight=transition_matrix[i, j])

    # Define positions for better layout
    pos = nx.spring_layout(G, seed=42)  # Spring layout for visualization

    # Draw the graph
    plt.figure(figsize=(10, 8))
    nx.draw(
        G, pos, with_labels=True, node_size=2000, node_color="lightblue", font_size=10, font_weight="bold"
    )

    # Draw edge labels with probabilities
    edge_labels = {
        (u, v): f"{d['weight']*100:.5f} %" for u, v, d in G.edges(data=True)
    }
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

    # Title and display the graph
    plt.title("Markov Transition Matrix Visualization", fontsize=14)
    plt.show()
