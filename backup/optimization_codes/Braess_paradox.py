# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 10:38:14 2024

@author: xzhou
"""


from pyomo.environ import *
from pyomo.opt import SolverFactory

# Define model
model = ConcreteModel()

# Define nodes and links for Braess's Paradox network structure
model.nodes = Set(initialize=['Start', 'A', 'B', 'End'])
model.links = Set(initialize=[('Start', 'A'), ('A', 'End'), ('Start', 'B'), ('B', 'End'), ('A', 'B')])

# Define parameters for free-flow time and capacity in minutes and vehicle units
free_flow_time = {
    ('Start', 'A'): 0.01,  # Represents T/100 for variable travel time based on flow
    ('A', 'End'): 45,      # Constant travel time link
    ('Start', 'B'): 45,    # Constant travel time link
    ('B', 'End'): 0.01,    # Represents T/100 for variable travel time based on flow
    ('A', 'B'): 0.001      # Very low travel time to demonstrate paradox effect
}

# Capacity is sufficiently high to simplify the effect of flow on time
capacity = {
    ('Start', 'A'): 4000,
    ('A', 'End'): 4000,
    ('Start', 'B'): 4000,
    ('B', 'End'): 4000,
    ('A', 'B'): 2000
}

# BPR parameters for calculating the cost as per flow and capacity constraints
alpha = 0.15
beta = 4

# Demand and node indicators for the flow balance
demand_vehicle = 4000  # Total vehicles from Start to End
origin_node = {'Start': 1, 'A': 0, 'B': 0, 'End': 0}
destination_node = {'Start': 0, 'A': 0, 'B': 0, 'End': 1}

# Define Pyomo variables
model.x = Var(model.links, domain=NonNegativeReals)  # Flow on each link

# User Equilibrium objective (minimizing individual cost)
def user_equilibrium_objective(model):
    return sum(
        free_flow_time[i, j] * (
            model.x[i, j] + (alpha * model.x[i, j] ** (beta + 1)) / capacity[i, j] ** beta
        )
        for i, j in model.links
    )

model.UE_cost = Objective(rule=user_equilibrium_objective, sense=minimize)

# System Optimal objective (minimizing total system cost)
def system_optimal_objective(model):
    return sum(
        free_flow_time[i, j] * model.x[i, j] * (1 + alpha * (model.x[i, j] / capacity[i, j]) ** beta)
        for i, j in model.links
    )

model.SO_cost = Objective(rule=system_optimal_objective, sense=minimize)

# Flow balance constraints to ensure conservation of vehicles at nodes
def flow_balance_rule(model, node):
    inflow = sum(model.x[i, j] for i, j in model.links if j == node)
    outflow = sum(model.x[i, j] for i, j in model.links if i == node)
    return inflow - outflow == (demand_vehicle if origin_node[node] else 0) - (demand_vehicle if destination_node[node] else 0)

model.flow_balance = Constraint(model.nodes, rule=flow_balance_rule)

# Define the solver and solve both objectives separately
solver = SolverFactory('ipopt')

# Solve User Equilibrium
model.SO_cost.deactivate()
result_UE = solver.solve(model, tee=True)
model.UE_cost.display()

# Solve System Optimal
model.SO_cost.activate()
model.UE_cost.deactivate()
result_SO = solver.solve(model, tee=True)
model.SO_cost.display()

# Results Summary
print("\nBraess's Paradox Model Results")
print("User Equilibrium Results (UE):")
for i, j in model.links:
    print(f"Flow on link ({i},{j}):", model.x[i, j].value)

print("\nSystem Optimal Results (SO):")
for i, j in model.links:
    print(f"Flow on link ({i},{j}):", model.x[i, j].value)
