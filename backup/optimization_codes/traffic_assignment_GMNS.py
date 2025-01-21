import pandas as pd
from pyomo.environ import *
from pyomo.opt import SolverFactory

# Load data from CSV files
nodes_df = pd.read_csv('node.csv')
links_df = pd.read_csv('link.csv')
demand_df = pd.read_csv('demand.csv')

# Prepare Sets and Parameters from the data

# Nodes set
nodes = set(nodes_df['node_id'])

# Links set and parameters
links = set(zip(links_df['from_node_id'], links_df['to_node_id']))
free_flow_time = {(row['from_node_id'], row['to_node_id']): row['length'] / row['free_speed'] for _, row in links_df.iterrows()}
capacity = {(row['from_node_id'], row['to_node_id']): row['capacity'] for _, row in links_df.iterrows()}
alpha = {(row['from_node_id'], row['to_node_id']): row['VDF_alpha'] for _, row in links_df.iterrows()}
beta = {(row['from_node_id'], row['to_node_id']): row['VDF_beta'] for _, row in links_df.iterrows()}

# Demand data for OD pairs
demand_data = {(row['o_zone_id'], row['d_zone_id']): row['volume'] for _, row in demand_df.iterrows() if row['volume'] > 0}

# Define the Pyomo model
model = ConcreteModel()

# Define Sets
model.nodes = Set(initialize=nodes)  # Nodes set
model.links = Set(initialize=links)  # Links set
model.commodities = Set(initialize=demand_data.keys())  # Set of (origin, destination) pairs

# Define Parameters
model.free_flow_time = Param(model.links, initialize=free_flow_time)
model.capacity = Param(model.links, initialize=capacity)
model.alpha = Param(model.links, initialize=alpha)
model.beta = Param(model.links, initialize=beta)
model.demand = Param(model.commodities, initialize=demand_data)

# Define Variables
# Flow on each link for each commodity
model.x = Var(model.commodities, model.links, domain=NonNegativeReals)

# Objective function: Minimize the total travel time cost across all commodities
def objective_rule(model):
    return sum(
        model.free_flow_time[i, j] * (
            sum(model.x[k, i, j] for k in model.commodities) + 
            (model.alpha[i, j] * sum(model.x[k, i, j] for k in model.commodities) ** (model.beta[i, j] + 1)) / 
            (model.capacity[i, j] ** model.beta[i, j] * (model.beta[i, j] + 1))
        )
        for (i, j) in model.links
    )
model.obj = Objective(rule=objective_rule, sense=minimize)

# Flow conservation constraints for each commodity
def flow_conservation_rule(model, k, node):
    o, d = k  # origin and destination for the current commodity
    inflow = sum(model.x[k, i, j] for (i, j) in model.links if j == node)
    outflow = sum(model.x[k, i, j] for (i, j) in model.links if i == node)
    
    if node == o:  # Origin node for this commodity
        return outflow - inflow == model.demand[k]
    elif node == d:  # Destination node for this commodity
        return inflow - outflow == model.demand[k]
    else:  # Intermediate nodes must have balanced flow
        return inflow == outflow

# Apply the rule to the constraint with proper indexing for commodities and nodes
model.flow_conservation = Constraint(model.commodities, model.nodes, rule=flow_conservation_rule)

# Solve the model
solver = SolverFactory('ipopt')  # You may use other solvers like 'glpk' or 'cbc'
solver.solve(model, tee=True)

# Display results
print("\nOptimal Flows on Links by Commodity:")
for k in model.commodities:
    for (i, j) in model.links:
        print(f"Flow of commodity {k} on link ({i}, {j}): {model.x[k, i, j].value:.4f} vehicles")

print(f"\nObjective value (Total Cost): {model.obj():.4f}")

# Compute and display travel times for each link using BPR function
def compute_travel_time(link_flow, link_capacity, free_flow_time, alpha, beta):
    return free_flow_time * (1 + alpha * (link_flow / link_capacity) ** beta)

print("\nLink Travel Times (using BPR function):")
for (i, j) in model.links:
    total_flow = sum(model.x[k, i, j].value for k in model.commodities)
    travel_time = compute_travel_time(total_flow, model.capacity[(i, j)], model.free_flow_time[(i, j)],
                                      model.alpha[(i, j)], model.beta[(i, j)])
    print(f"Travel time on link ({i}, {j}): {travel_time:.4f} minutes")
