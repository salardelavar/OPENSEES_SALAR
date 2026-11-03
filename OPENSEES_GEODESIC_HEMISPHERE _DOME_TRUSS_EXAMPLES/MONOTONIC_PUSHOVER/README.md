# MONOTONIC PUSHOVER ANALYSIS OF GEODESIC HEMISPHERE DOME TRUSS WITH GEOMETRIC AND MATERIAL NONLINEARITY IN OPENSEES

This OpenSeesPy script performs a monotonic displacement-controlled pushover analysis of a geodesic hemisphere dome truss that accounts for both geometric and material nonlinearity. It first generates the geometry by creating a frequency-subdivided icosahedron, extracting the upper hemisphere nodes and truss edges, then builds a 3D OpenSees model in which base nodes are fixed, steel members are assigned either elastic or hysteretic inelastic material, and either ordinary Truss or CorotTruss (large-displacement) elements are used. A linearly varying vertical (or horizontal) reference load is applied to all non-base nodes, after which a DisplacementControl integrator drives the apex node to a prescribed target displacement while recording the corresponding base reaction; finally the script plots the original versus deformed geometry together with the resulting pushover curve.

![alt text](geodesic_4v_pushover_DC)

# Multilinear Fitting of Monotonic Pushover Curve
![alt text](Multilinear_Fitting_of_Monotonic_Pushover_Curve)
