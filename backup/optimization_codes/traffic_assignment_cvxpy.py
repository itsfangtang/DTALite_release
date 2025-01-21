# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 21:51:01 2024

@author: xzhou
"""

# -*- coding: utf-8 -*-
"""
Traffic Assignment Problem User Equilibrium using CVXPY

Author: xzhou
Created on Sat Nov  2 19:41:58 2024
"""

import cvxpy as cp
import numpy as np

# Define the sets
nodes = [1, 2, 3, 4]
links = [(1, 3), (3, 2), (1, 4), (4, 2)]

# Define parameters
free_flow_time = {
    (1, 3): 10,
    (3, 2): 10,
    (1, 4): 15,
    (4, 2): 15,
}

capacity = {
    (1, 3): 4000,
    (3, 2): 4000,
    (1, 4): 3000,
    (4, 2): 4000,
}

# BPR parameters
alpha = 0.15
beta = 2

# Demand and node indicators
demand_vehicle = 1000
origin_node = {1: 1, 2: 0, 3: 0, 4: 0}
destination_node = {1: 0, 2: 1, 3: 0, 4: 0}

# Intermediate node indicator
intermediate_node = {i: (1 - origin_node[i]) * (1 - destination_node[i]) for i in nodes}

# Define CVXPY variables
x = {link: cp.Variable(nonneg=True) for link in links}  # Flow on each link

# Define the objective function
objective_terms = []
for (i, j) in links:
    term = free_flow_time[(i, j)] * (x[(i, j)] + 
           (alpha * cp.power(x[(i, j)], beta + 1)) / (capacity[(i, j)]**beta * (beta + 1)))
    objective_terms.append(term)

objective = cp.Minimize(cp.sum(objective_terms))

# Define constraints
constraints = []

# Flow conservation constraints with detailed prints for debugging
for node in nodes:
    inflow = []
    outflow = []
    for (i, j) in links:
        if j == node:
            inflow.append(x[(i, j)])
        if i == node:
            outflow.append(x[(i, j)])
    
    # Print the intermediate node indicator for debugging
    print(f"Node {node}:")
    print(f"  Intermediate Node Indicator: {intermediate_node[node]}")
    print(f"  Inflow: {[str(inflow_expr) for inflow_expr in inflow]}")
    print(f"  Outflow: {[str(outflow_expr) for outflow_expr in outflow]}")

    # Check node type and append corresponding constraint
    if origin_node[node] == 1:
        # Origin node: Outflow - Inflow = Demand
        constraint = cp.sum(outflow) - cp.sum(inflow) == demand_vehicle
        constraints.append(constraint)
        print(f"  Constraint (Origin): {constraint}")
    elif destination_node[node] == 1:
        # Destination node: Inflow - Outflow = Demand
        constraint = cp.sum(inflow) - cp.sum(outflow) == demand_vehicle
        constraints.append(constraint)
        print(f"  Constraint (Destination): {constraint}")
    elif intermediate_node[node] == 1:
        # Intermediate nodes: Inflow == Outflow
        constraint = cp.sum(inflow) == cp.sum(outflow)
        constraints.append(constraint)
        print(f"  Constraint (Intermediate): {constraint}")
    else:
        print("  No constraint added for this node.")

  

# Form and solve the problem
prob = cp.Problem(objective, constraints)
prob.solve(solver=cp.SCS, verbose=True)  # You can choose other solvers like 'ECOS', 'OSQP'

# Display results
print("\nOptimal Flows on Links:")
for (i, j) in links:
    flow_value = x[(i, j)].value
    if flow_value is not None:
        print(f"Flow on link ({i}, {j}): {flow_value:.4f} vehicles")
    else:
        print(f"Flow on link ({i}, {j}): No solution found")


print(f"\nObjective value (Total Cost): {prob.value:.4f}")



