from pyomo.environ import *
from pyomo.opt import SolverFactory

# Define model
model = ConcreteModel()

# Define the sets
model.nodes = Set(initialize=[1, 2, 3, 4])
model.links = Set(initialize=[(1, 3), (3, 2), (1, 4), (4, 2)])

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
beta = 4

# Demand and node indicators
demand_vehicle = 8000
origin_node = {1: 1, 2: 0, 3: 0, 4: 0}
destination_node = {1: 0, 2: 1, 3: 0, 4: 0}

# Intermediate node indicator
intermediate_node = {i: (1 - origin_node[i]) * (1 - destination_node[i]) for i in model.nodes}

# Define Pyomo variables
model.x = Var(model.links, domain=NonNegativeReals)  # Flow on each link

# Objective function components
def objective_rule(model):
    return sum(
        free_flow_time[i, j] * (
            model.x[i, j] + (alpha * model.x[i, j] ** (beta + 1)) / (capacity[i, j] ** beta * (beta + 1))
        )
        for (i, j) in model.links
    )

model.obj = Objective(rule=objective_rule, sense=minimize)

# Flow conservation constraints
def flow_conservation_rule(model, node):
    inflow = sum(model.x[i, j] for (i, j) in model.links if j == node)
    outflow = sum(model.x[i, j] for (i, j) in model.links if i == node)
    
    if origin_node[node] == 1:
        return outflow - inflow == demand_vehicle
    elif destination_node[node] == 1:
        return inflow - outflow == demand_vehicle
    elif intermediate_node[node] == 1:
        return inflow == outflow
    else:
        return Constraint.Skip

model.flow_conservation = Constraint(model.nodes, rule=flow_conservation_rule)

# Solve the problem
solver = SolverFactory('ipopt')  # You can use other solvers like 'glpk', 'cbc', 'gurobi', etc.
solver.solve(model, tee=True)

# Display results
print("\nOptimal Flows on Links:")
for (i, j) in model.links:
    print(f"Flow on link ({i}, {j}): {model.x[i, j].value:.4f} vehicles")

print(f"\nObjective value (Total Cost): {model.obj():.4f}")


# Calculate and display path flows
# Path 1: 1 -> 3 -> 2
path1_flow = min(model.x[1, 3].value, model.x[3, 2].value)
# Path 2: 1 -> 4 -> 2
path2_flow = min(model.x[1, 4].value, model.x[4, 2].value)

print("\nPath Flows:")
print(f"Flow on Path 1 (1->3->2): {path1_flow:.4f} vehicles")
print(f"Flow on Path 2 (1->4->2): {path2_flow:.4f} vehicles")

# Compute link travel times using BPR function
def compute_travel_time(link_flow, link_capacity, t0):
    return t0 * (1 + alpha * (link_flow / link_capacity) ** beta)

travel_time_13 = compute_travel_time(model.x[1, 3].value, capacity[(1, 3)], free_flow_time[(1, 3)])
travel_time_32 = compute_travel_time(model.x[3, 2].value, capacity[(3, 2)], free_flow_time[(3, 2)])
travel_time_14 = compute_travel_time(model.x[1, 4].value, capacity[(1, 4)], free_flow_time[(1, 4)])
travel_time_42 = compute_travel_time(model.x[4, 2].value, capacity[(4, 2)], free_flow_time[(4, 2)])

# Compute path travel times
path1_travel_time = travel_time_13 + travel_time_32
path2_travel_time = travel_time_14 + travel_time_42

print("\nPath Travel Times:")
print(f"Travel Time on Path 1 (1->3->2): {path1_travel_time:.4f} minutes")
print(f"Travel Time on Path 2 (1->4->2): {path2_travel_time:.4f} minutes")

# Determine the minimum travel time
min_travel_time = min(path1_travel_time, path2_travel_time)
print(f"\nMinimum Travel Time: {min_travel_time:.4f} minutes")

# Compute User Equilibrium Gaps
# In UE, all used paths should have travel time equal to the minimum travel time
gap_path1 = path1_travel_time - min_travel_time
gap_path2 = path2_travel_time - min_travel_time

print("\nUser Equilibrium Gaps:")
print(f"Gap for Path 1 (1->3->2): {gap_path1:.4f} minutes")
print(f"Gap for Path 2 (1->4->2): {gap_path2:.4f} minutes")


# Weighted Gaps: Gap * Flow
weighted_gap_path1 = path1_flow * gap_path1
weighted_gap_path2 = path2_flow * gap_path2

# Total Weighted Gap
total_weighted_gap = weighted_gap_path1 + weighted_gap_path2

print("\nUser Equilibrium Gaps (Weighted by Flow):")
print(f"Weighted Gap for Path 1 (1->3->2): {weighted_gap_path1:.4f} vehicle-minutes")
print(f"Weighted Gap for Path 2 (1->4->2): {weighted_gap_path2:.4f} vehicle-minutes")
print(f"Total Weighted Gap: {total_weighted_gap:.4f} vehicle-minutes")

# Optionally, display the objective value
print(f"\nObjective value (Total Cost): {model.obj():.4f}")
