# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 12:55:53 2024

@author: xzhou
"""

from pyomo.environ import *
from pyomo.opt import SolverFactory
import pandas as pd

# Define the Pyomo model
model = ConcreteModel()

# Sets
model.i = Set(initialize=[1, 2, 3])  # Nodes
model.t = Set(initialize=[0, 1, 2, 3, 4, 5])  # Time steps
model.a = Set(initialize=[1, 2, 3])  # Agents
model.k = Set(initialize=[1, 2])  # Additional nodes (for constraints)

# Parameters
c = {(1, 2, 0, 1): 1, (1, 2, 2, 3): 1, (2, 3, 1, 3): 2, (1, 3, 1, 4): 3}  # Sample travel costs
net = {(1, 2, 0, 1): 1, (1, 2, 2, 3): 1, (2, 3, 1, 3): 1, (1, 3, 1, 4): 1}
Cap = {(1, 2, 0, 1): 1, (1, 2, 2, 3): 2, (2, 3, 1, 3): 1, (1, 3, 1, 4): 1}
ccost = {(1, 2): 0, (2, 3): 0, (1, 3): 10}  # Construction costs
Budget = 15
original_least_path_time = 3
M = 10
epsilon = {1: 0, 2: 2, 3: 2}
r = {(1, 1): 1, (2, 2): 1, (1, 2): 1, (2, 3): 1, (1, 3): 0}  # Existing link factor

# Intermediate node indicator, origin and destination based on time
origin = {(1, 1, 0): 1, (2, 1, 0): 1, (3, 1, 0): 1}
destination = {(1, 3, 5): 1, (2, 3, 5): 1, (3, 3, 5): 1}
node_index = {(1, 1): 1, (2, 1): 1, (2, 3): 1, (3, 4): 1}
intermediate = {(a, i, t): (1 - origin.get((a, i, t), 0)) * (1 - destination.get((a, i, t), 0)) for a in model.a for i in model.i for t in model.t}

# Initialize parameters as Pyomo Parameters
model.c = Param(model.i, model.i, model.t, model.t, initialize=c, default=0)
model.net = Param(model.i, model.i, model.t, model.t, initialize=net, default=0)
model.Cap = Param(model.i, model.i, model.t, model.t, initialize=Cap, default=0)
model.ccost = Param(model.i, model.i, initialize=ccost, default=0)
model.Budget = Param(initialize=Budget)
model.original_least_path_time = Param(initialize=original_least_path_time)
model.M = Param(initialize=M)
model.epsilon = Param(model.a, initialize=epsilon)
model.r = Param(model.i, model.i, initialize=r, default=0)
model.origin = Param(model.a, model.i, model.t, initialize=origin, default=0)
model.destination = Param(model.a, model.i, model.t, initialize=destination, default=0)
model.node_index = Param(model.i, model.t, initialize=node_index, default=0)

# Variables
model.z = Var()
model.z_1 = Var()
model.z_2 = Var()
model.f = Var(model.i, model.i, model.t, model.t, domain=NonNegativeReals)
model.x = Var(model.a, model.i, model.i, model.t, model.t, domain=NonNegativeReals)
model.y = Var(model.i, model.i, domain=Binary)

# Fix binary variables for existing links
for (i, j), value in r.items():
    if value > 0:
        model.y[i, j].fix(1)

# Objective function
def so_obj_rule(model):
    return model.z == sum(model.f[i, j, t, s] * model.c[i, j, t, s] for i in model.i for j in model.i for t in model.t for s in model.t if model.net[i, j, t, s] > 0)
model.so_obj = Constraint(rule=so_obj_rule)

# Constraints for node flow conservation and capacity
def flow_conservation_rule(model, a, i, t):
    if model.origin[a, i, t] == 1:
        return sum(model.x[a, i, j, t, s] for j in model.i for s in model.t if model.net[i, j, t, s] > 0) == model.origin[a, i, t]
    elif model.destination[a, i, t] == 1:
        return sum(model.x[a, j, i, s, t] for j in model.i for s in model.t if model.net[j, i, s, t] > 0) == model.destination[a, i, t]
    elif intermediate.get((a, i, t), 0) == 1:
        return (sum(model.x[a, i, j, t, s] for j in model.i for s in model.t if model.net[i, j, t, s] > 0) - 
                sum(model.x[a, j, i, s, t] for j in model.i for s in model.t if model.net[j, i, s, t] > 0)) == 0
    return Constraint.Skip
model.flow_conservation = Constraint(model.a, model.i, model.t, rule=flow_conservation_rule)

# Budget constraint
def budget_constraint_rule(model):
    return sum(model.ccost[i, j] * model.y[i, j] for i in model.i for j in model.i) <= model.Budget
model.budget_constraint = Constraint(rule=budget_constraint_rule)

# Link capacity constraints
def capacity_constraint_rule(model, i, j, t, s):
    if model.net[i, j, t, s] > 0:
        return model.Cap[i, j, t, s] * model.y[i, j] - model.f[i, j, t, s] >= 0
    return Constraint.Skip
model.capacity_constraint = Constraint(model.i, model.i, model.t, model.t, rule=capacity_constraint_rule)

# Solve the model
solver = SolverFactory('scip')  # or any other solver like 'glpk' or 'cbc'
solver.solve(model, tee=True)

# Display results
for v in model.component_objects(Var, active=True):
    print(f"{v.name} = {v.value}")
