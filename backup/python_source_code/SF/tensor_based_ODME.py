import pandas as pd
import numpy as np
import tensorflow as tf
import os

# Function to process network data
def process_network_data(node_file, link_file, route_assignment_file):
    """
    Processes network data to generate matrices A_PL and B_OD_P.
    
    Args:
    - node_file: Path to the node.csv file (not used in this function).
    - link_file: Path to the link.csv file.
    - route_assignment_file: Path to the route_assignment.csv file.

    Returns:
    - A_PL: Path-to-Link incidence matrix.
    - B_OD_P: OD-to-Path incidence matrix.
    """
    link_df = pd.read_csv(link_file)
    route_assignment_df = pd.read_csv(route_assignment_file)

    # Create Path-to-Link Matrix (A_PL)
    num_paths = len(route_assignment_df)
    num_links = len(link_df)
    A_PL = np.zeros((num_paths, num_links), dtype=np.float32)
    for i, row in route_assignment_df.iterrows():
        link_ids = [int(link_id) for link_id in str(row['link_sequence']).split(';') if link_id.strip()]
        for link_id in link_ids:
            A_PL[i, link_id - 1] = 1  # Use 1-based indexing for link_id

    # Create OD-to-Path Matrix (B_OD_P)
    od_pairs = route_assignment_df[['o_zone_id', 'd_zone_id']].drop_duplicates().values
    num_od_pairs = len(od_pairs)
    B_OD_P = np.zeros((num_od_pairs, num_paths), dtype=np.float32)
    for i, (origin, destination) in enumerate(od_pairs):
        matching_paths = route_assignment_df[
            (route_assignment_df['o_zone_id'] == origin) & (route_assignment_df['d_zone_id'] == destination)
        ].index
        for path_id in matching_paths:
            B_OD_P[i, path_id] = 1

    return tf.constant(A_PL, dtype=tf.float32), tf.constant(B_OD_P, dtype=tf.float32)

# Function to read link data
def read_link_file(link_file):
    """
    Reads the link.csv file to extract observed link flows and attributes.
    
    Args:
    - link_file: Path to the link.csv file.

    Returns:
    - f_L_obs: Observed link flows (ref_volume column).
    - C_L: Link capacities (capacity column).
    - T_L_0: Free-flow travel times (free_speed column).
    """
    link_df = pd.read_csv(link_file)
    if 'ref_volume' not in link_df.columns:
        raise ValueError("The 'ref_volume' column is missing in link.csv.")
    
    f_L_obs = link_df['ref_volume'].values.astype(np.float32)
    C_L = link_df['capacity'].values.astype(np.float32)
    T_L_0 = link_df['free_speed'].values.astype(np.float32)
    
    return tf.constant(f_L_obs), tf.constant(C_L), tf.constant(T_L_0)

# Function to read demand data
def read_demand_file(demand_file, route_assignment):
    """
    Reads the demand.csv file to extract observed OD travel times or volumes.
    
    Args:
    - demand_file: Path to the demand.csv file.

    Returns:
    - T_OD_obs: Observed OD travel times or volumes.
    """
    demand_df = pd.read_csv(demand_file)
    route_assignment_df = pd.read_csv(route_assignment_file)
    init_od_df = route_assignment_df[['o_zone_id', 'd_zone_id']].drop_duplicates()

    if init_od_df.shape[0] != demand_df.shape[0]:
        print(f"WARNING: The length of target OD flows is {demand_df.shape[0]}, "
              f"but the length of initial OD flows is {init_od_df.shape[0]}. ")
        print()
        print("Imputing the od flow target data to remove unmatched OD pairs...")
        print(f" - Before removing, the unmatched target od volume is {demand_df['volume'].sum()}", )
        demand_df = \
            demand_df.merge(init_od_df[["o_zone_id", "d_zone_id"]], on=["o_zone_id", "d_zone_id"], how="inner")
        print(f" - After processing, the target od volume is {demand_df['volume'].sum()}")
    else:
        print("The lengths of target and initial OD flows match, no action required.")

    if 'volume' not in demand_df.columns:
        raise ValueError("The 'volume' column is missing in demand.csv.")
    
    T_OD_obs = demand_df['volume'].values.astype(np.float32)
    return tf.constant(T_OD_obs)

# OD Demand Estimation Function
def od_demand_estimation_with_observations(A_PL, B_OD_P, C_L, T_L_0, f_L_obs, T_OD_obs,
                                           alpha=0.15, beta=4.0, num_steps=500, lr=0.1,
                                           lambda_1=1.0, lambda_2=1.0):
    """
    Estimate OD demand using observed link flows (f_L_obs) and OD travel times (T_OD_obs).
    
    Args:
    - A_PL: Path-to-Link incidence matrix.
    - B_OD_P: OD-to-Path incidence matrix.
    - C_L: Link capacities.
    - T_L_0: Free-flow travel times for links.
    - f_L_obs: Observed link flows.
    - T_OD_obs: Observed OD travel times or volumes.
    - alpha: BPR function parameter (default: 0.15).
    - beta: BPR function parameter (default: 4.0).
    - num_steps: Number of optimization steps (default: 500).
    - lr: Learning rate for the optimizer (default: 0.1).
    - lambda_1, lambda_2: Weights for link flow and OD travel time deviations.

    Returns:
    - Estimated OD demands (f_OD).
    """
    num_od_pairs = B_OD_P.shape[0]

    # Initialize OD Demand (f_OD)
    f_OD = tf.Variable(tf.ones(num_od_pairs, dtype=tf.float32) * 100.0)  # Initial guess

    # Objective Function
    def objective_function():
        # Step 1: Compute Path Flows (f_P)
        f_P = tf.matmul(tf.transpose(B_OD_P), tf.expand_dims(f_OD, axis=1))[:, 0]

        # Step 2: Compute Link Flows (f_L)
        f_L = tf.matmul(tf.transpose(A_PL), tf.expand_dims(f_P, axis=1))[:, 0]

        # Step 3: Compute Link Travel Times (T_L)
        T_L = T_L_0 + alpha * (f_L / C_L) ** beta

        # Step 4: Compute Path Travel Times (T_P)
        T_P = tf.matmul(A_PL, tf.expand_dims(T_L, axis=1))[:, 0]

        # Step 5: Compute OD Travel Times (T_OD)
        T_OD = tf.matmul(B_OD_P, tf.expand_dims(T_P, axis=1))[:, 0]

        # Deviations from Observations
        loss_link_flows = lambda_1 * tf.reduce_sum((f_L_obs - f_L) ** 2)
        loss_od_travel_times = lambda_2 * tf.reduce_sum((T_OD_obs - T_OD) ** 2)
        
        # Total Loss
        total_loss = loss_link_flows + loss_od_travel_times
        return total_loss

    # Optimizer
    optimizer = tf.keras.optimizers.Adam(learning_rate=lr)

    # Training Loop
    for step in range(num_steps):
        with tf.GradientTape() as tape:
            loss = objective_function()

        # Compute and apply gradients
        gradients = tape.gradient(loss, [f_OD])
        optimizer.apply_gradients(zip(gradients, [f_OD]))

        # Enforce Non-Negativity
        f_OD.assign(tf.maximum(f_OD, 0.0))

        # Logging every 100 steps
        if step % 100 == 0 or step == num_steps - 1:
            print(f"Step {step}, Loss: {loss.numpy():.2f}")
            print(f"Estimated OD Demands (f_OD): {f_OD.numpy()}")

    return f_OD.numpy()

# Main Execution
# Define the file paths to access CSV files
current_path = os.getcwd()
upper_path = os.path.dirname(current_path)

# File Paths
node_file = os.path.join(upper_path, "node.csv") #"node.csv"  # Path to node file (not used in this function)
link_file = os.path.join(upper_path, "link.csv") #"link.csv"  # Path to link file
route_assignment_file = os.path.join(upper_path, "route_assignment.csv")# "route_assignment.csv"  # Path to route assignment file
demand_file = os.path.join(upper_path, "demand.csv") # "demand.csv"  # Path to demand file

# Process network data
A_PL, B_OD_P = process_network_data(node_file, link_file, route_assignment_file)

# Read link and demand data
f_L_obs, C_L, T_L_0 = read_link_file(link_file)
T_OD_obs = read_demand_file(demand_file, route_assignment_file)

# Run OD Demand Estimation
estimated_f_OD = od_demand_estimation_with_observations(A_PL, B_OD_P, C_L, T_L_0, f_L_obs, T_OD_obs)

# Print Final Estimated OD Demands
print("\nFinal Estimated OD Demands:", estimated_f_OD)
