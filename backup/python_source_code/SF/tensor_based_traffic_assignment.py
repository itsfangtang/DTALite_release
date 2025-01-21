# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 16:03:10 2024

@author: xzhou
"""

import tensorflow as tf
import pandas as pd
import numpy as np

# Function to process network data
def process_network_data(node_file, link_file, route_assignment_file):
    # Load the input files
    node_df = pd.read_csv(node_file)
    link_df = pd.read_csv(link_file)
    route_assignment_df = pd.read_csv(route_assignment_file)

    # Create Path-to-Link Matrix (A_PL)
    num_paths = len(route_assignment_df)
    num_links = len(link_df)
    A_PL = np.zeros((num_paths, num_links), dtype=np.float32)
    for i, row in route_assignment_df.iterrows():
        link_ids = [int(link_id) for link_id in str(row['link_ids']).split(';') if link_id.strip()]
        for link_id in link_ids:
            A_PL[i, link_id - 1] = 1  # Use 1-based indexing for link_id

    # Create OD-to-Path Matrix (B_OD_P)
    od_pairs = route_assignment_df[['origin', 'destination']].drop_duplicates().values
    num_od_pairs = len(od_pairs)
    B_OD_P = np.zeros((num_od_pairs, num_paths), dtype=np.float32)
    for i, (origin, destination) in enumerate(od_pairs):
        matching_paths = route_assignment_df[
            (route_assignment_df['origin'] == origin) & (route_assignment_df['destination'] == destination)
        ].index
        for path_id in matching_paths:
            B_OD_P[i, path_id] = 1

    # Extract OD demands
    od_demands = route_assignment_df.groupby(['origin', 'destination'])['total_travel_time'].sum().values

    # Extract Link Data: Capacities and Free-Flow Times
    C_L = link_df['capacity'].values.astype(np.float32)  # Link capacities
    T_L_0 = link_df['free_speed'].values.astype(np.float32)  # Free-flow travel times

    # Convert matrices to TensorFlow tensors
    A_PL_tensor = tf.constant(A_PL, dtype=tf.float32)
    B_OD_P_tensor = tf.constant(B_OD_P, dtype=tf.float32)
    od_demands_tensor = tf.constant(od_demands, dtype=tf.float32)
    C_L_tensor = tf.constant(C_L, dtype=tf.float32)
    T_L_0_tensor = tf.constant(T_L_0, dtype=tf.float32)

    return A_PL_tensor, B_OD_P_tensor, od_demands_tensor, C_L_tensor, T_L_0_tensor

# Traffic Assignment Function
def traffic_assignment(A_PL, B_OD_P, f_OD, C_L, T_L_0, alpha=0.15, beta=4.0, num_steps=500, lr=0.1):
    # Initialize Path Flows (f_P) evenly across paths
    num_paths = A_PL.shape[0]
    f_P = tf.Variable(tf.zeros(num_paths, dtype=tf.float32))

    def initialize_path_flows():
        path_counts = tf.reduce_sum(B_OD_P, axis=1)  # Count of paths for each OD pair
        initial_flows = tf.matmul(tf.transpose(B_OD_P), tf.expand_dims(f_OD / path_counts, axis=1))[:, 0]
        f_P.assign(initial_flows)

    initialize_path_flows()

    # Objective Function
    def objective_function():
        # Compute Link Flows (f_L)
        f_L = tf.matmul(tf.transpose(A_PL), tf.expand_dims(f_P, axis=1))[:, 0]

        # Compute Link Travel Times (T_L)
        T_L = T_L_0 + alpha * (f_L / C_L) ** beta

        # Compute Path Travel Times (T_P)
        T_P = tf.matmul(A_PL, tf.expand_dims(T_L, axis=1))[:, 0]

        # Compute OD Travel Times (T_OD)
        T_OD = tf.matmul(B_OD_P, tf.expand_dims(T_P, axis=1))[:, 0]

        # Total Travel Time
        total_travel_time = tf.reduce_sum(f_OD * T_OD)
        return total_travel_time

    # Constraint Enforcement for OD Flows
    def enforce_od_flow_constraints():
        current_f_OD = tf.matmul(B_OD_P, tf.expand_dims(f_P, axis=1))[:, 0]
        correction = current_f_OD - f_OD
        f_P.assign_sub(tf.matmul(tf.transpose(B_OD_P), tf.expand_dims(correction, axis=1))[:, 0])

    # Optimizer
    optimizer = tf.keras.optimizers.Adam(learning_rate=lr)

    # Training Loop
    for step in range(num_steps):
        with tf.GradientTape() as tape:
            loss = objective_function()

        # Compute and apply gradients
        gradients = tape.gradient(loss, [f_P])
        optimizer.apply_gradients(zip(gradients, [f_P]))

        # Enforce Non-Negativity
        f_P.assign(tf.maximum(f_P, 0.0))

        # Enforce OD Flow Constraints
        enforce_od_flow_constraints()

        # Logging every 100 steps
        if step % 100 == 0:
            print(f"Step {step}, Loss: {loss.numpy():.2f}")

    # Final Results
    f_L = tf.matmul(tf.transpose(A_PL), tf.expand_dims(f_P, axis=1))[:, 0]
    T_L = T_L_0 + alpha * (f_L / C_L) ** beta
    T_P = tf.matmul(A_PL, tf.expand_dims(T_L, axis=1))[:, 0]
    T_OD = tf.matmul(B_OD_P, tf.expand_dims(T_P, axis=1))[:, 0]

    return f_P.numpy(), f_L.numpy(), T_L.numpy(), T_P.numpy(), T_OD.numpy()

# Load the files (replace paths with your local file paths)
node_file = "node.csv"
link_file = "link.csv"
route_assignment_file = "route_assignment.csv"

# Process the network data
A_PL, B_OD_P, f_OD, C_L, T_L_0 = process_network_data(node_file, link_file, route_assignment_file)

# Run the traffic assignment
f_P, f_L, T_L, T_P, T_OD = traffic_assignment(A_PL, B_OD_P, f_OD, C_L, T_L_0)

# Print the final results
print("\nFinal Path Flows (f_P):", f_P)
print("Final Link Flows (f_L):", f_L)
print("Final Link Travel Times (T_L):", T_L)
print("Final Path Travel Times (T_P):", T_P)
print("Final OD Travel Times (T_OD):", T_OD)
