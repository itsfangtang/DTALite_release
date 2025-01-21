# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:43:15 2024

@author: xzhou
"""

import csv
import os
import math
import pandas as pd
import osm2gmns as og

from collections import defaultdict

import numpy as np
import time
from datetime import datetime

MAX_MODE_TYPES = 10
# Global Variables
number_of_zones = 0
number_of_modes = 1
number_of_nodes = 0
number_of_links = 0
max_routes = 10;
first_thru_node = 0

AssignIterations = 20
demand_period_starting_hours = 7
demand_period_ending_hours = 8
g_tap_log_file = 0
g_base_demand_mode = 1
g_ODME_mode = 0
g_ODME_obs_VMT = -1
g_System_VMT = 0

g_ODME_link_volume_penalty = 0.01  # Relative weight on volume, converts deviation of link volume to travel time
g_ODME_VMT_penalty = 0.01

# Global dictionaries
g_map_external_node_id_2_node_seq_no = {}  # Maps external node ID to node sequence number
g_map_external_node_id_2_zone_id = {}  # Maps external node ID to zone_id

g_map_node_seq_no_2_external_node_id = {}  # Maps node sequence number to external node ID

g_link_vector = []  # A list to hold LinkRecord objects

# Replace `int* FirstLinkFrom;` and `int* LastLinkFrom;` with Python lists
FirstLinkFrom = []  # A list to store the first link from each node
LastLinkFrom = []   # A list to store the last link from each node

# Replace `sorted_list* LinksTo;` with a dictionary of lists for better handling
LinksTo = {}  # A dictionary where keys are node IDs and values are lists of g_link_vector directed to the node

# Global 5D list for storing link sequences
g_link_indices = []

# Placeholder definitions

total_o_flow = None       # Placeholder for a 1D array
g_zone_outbound_link_size = None  # Placeholder for a 1D integer array
md_od_flow = None         # Placeholder for a 3D array
md_route_cost = None      # Placeholder for a 3D array

def osm2gmns_network():
    
    input_file = 'map.osm'
    net = og.getNetFromFile(input_file, link_types=('motorway','trunk','primary','secondary')) 
    
    og.consolidateComplexIntersections(net, auto_identify=True)
    og.fillLinkAttributesWithDefaultValues(net, default_lanes=True, default_speed=True, default_capacity=True)
    og.generateNodeActivityInfo(net)
    og.outputNetToCSV(net)

# Function to calculate Euclidean distance between two nodes
def calculate_distance(node1, node2):
    return math.sqrt((node1["x_coord"] - node2["x_coord"])**2 + (node1["y_coord"] - node2["y_coord"])**2)

# Read and process the `node.csv` file
def trip_generation():
    file_path  = 'node.csv'
    node_list = []  # List to store all nodes
    activity_nodes = []  # List to store nodes with a positive zone_id
    production_totals = {}  # Dictionary to store production totals per node
    number_of_zones = 0  # Track the maximum zone_id

    # Step 1: Read the node file
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        for row in csv_reader:
            # Extract node_id, zone_id, and coordinates
            node_id = int(row.get("node_id", 0))
            zone_id_value = row.get("zone_id", "")
            zone_id = int(zone_id_value) if zone_id_value.isdigit() else 0

            # Read x_coord and y_coord
            x_coord = float(row.get("x_coord", 0))
            y_coord = float(row.get("y_coord", 0))

            # Add node to node list
            node = {
                "node_id": node_id,
                "zone_id": zone_id,
                "x_coord": x_coord,
                "y_coord": y_coord,
                "production": 1,  # Optional initial production
                "attraction": 1  # Optional initial attraction
            }
            node_list.append(node)
            
            # Add to activity nodes if zone_id > 0
            if zone_id > 0:
                activity_nodes.append(node)
                # Update the maximum zone_id
                if zone_id > number_of_zones:
                    number_of_zones = zone_id

            # Initialize production totals for each node as 1 (including itself)
            production_totals[node_id] = 0

    # Step 2: Map nodes to the nearest activity node
    for node in node_list:
        if node["zone_id"] > 0:
            continue  # Skip activity nodes
        
        min_distance = float("inf")
        nearest_activity_node = None

        for activity_node in activity_nodes:
            distance = calculate_distance(node, activity_node)
            if distance < min_distance:
                min_distance = distance
                nearest_activity_node = activity_node
        
        # Update the production total of the nearest activity node
        if nearest_activity_node:
            production_totals[nearest_activity_node["node_id"]] += 1

    # Step 3: Update the production column for activity nodes
    updated_nodes = []
    for node in node_list:
        if node["node_id"] in production_totals:
            node["production"] = production_totals[node["node_id"]]
        updated_nodes.append(node)

    # Step 4: Rewrite the `node.csv` file with updated columns
    with open(file_path, 'w', newline='') as csv_file:
        fieldnames = ["node_id", "zone_id", "x_coord", "y_coord", "production", "attraction"]
        csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        csv_writer.writeheader()
        for node in updated_nodes:
            csv_writer.writerow({
                "node_id": node["node_id"],
                "zone_id": node["zone_id"],
                "x_coord": node["x_coord"],
                "y_coord": node["y_coord"],
                "production": node["production"]*100,
                "attraction": node["production"]*100  # for simplicity
            })

def perform_simple_logit_model_for_trip_distribution(beta=0.1):
    """
    Reads `node.csv` and `accessibility_matrix.csv`, performs a simple logit-based destination choice model,
    and writes the output to `demand.csv`.

    :param beta: Coefficient for the utility function using accessibility
    """
    import pandas as pd
    import math

    node_file = 'node.csv'
    accessibility_matrix_file = 'accessibility_matrix.csv'
    output_demand_file = 'demand.csv'

    # Step 1: Read node.csv
    node_df = pd.read_csv(node_file)

    # Filter for zones with zone_id > 0 (activity zones)
    zones_df = node_df[node_df['zone_id'] > 0]

    # Create dictionaries for production and attraction values by o_zone_id and d_zone_id
    production = zones_df.set_index('zone_id')['production'].to_dict()
    attraction = zones_df.set_index('zone_id')['attraction'].to_dict()

    # Step 2: Read accessibility_matrix.csv
    accessibility_df = pd.read_csv(accessibility_matrix_file)

    # Create an empty list to store demand data
    demand_data = []

    # Step 3: Perform destination choice model
    for o_zone_id in production:
        # Filter rows in the accessibility matrix for the current origin zone
        o_accessibility = accessibility_df[accessibility_df['o_zone_id'] == o_zone_id]

        # Calculate total utility for normalization
        total_utility = 0
        utility_dict = {}
        for _, row in o_accessibility.iterrows():
            d_zone_id = row['d_zone_id']
            accessibility = row['cost']

            # Skip if destination is not in the attraction dictionary
            if d_zone_id not in attraction:
                continue

            # Compute utility
            utility = math.exp(-beta * accessibility)
            utility_dict[d_zone_id] = utility
            total_utility += utility

        # Compute demand for each destination zone
        for d_zone_id, utility in utility_dict.items():
            if total_utility > 0:  # Avoid division by zero
                prob = utility / total_utility  # Logit model
                Tij = production[o_zone_id] * prob

                # Append to demand data
                demand_data.append({
                    'o_zone_id': o_zone_id,
                    'd_zone_id': d_zone_id,
                    'volume': Tij
                })

    # Step 4: Write demand.csv
    demand_df = pd.DataFrame(demand_data)
    # Ensure o_zone_id and d_zone_id columns are integers
    demand_df['o_zone_id'] = demand_df['o_zone_id'].astype(int)
    demand_df['d_zone_id'] = demand_df['d_zone_id'].astype(int)
    
    demand_df.to_csv(output_demand_file, index=False)

    print(f"Demand matrix has been successfully written to {output_demand_file}.")
      
 

def perform_gravity_model(beta=0.1, model_type="singly"):
    node_file = 'node.csv'
    accessibility_matrix_file =  "accessibility_matrix.csv"
    output_demand_file = 'demand.csv'
    """
    Reads `node.csv` and `accessibility_matrix.csv`, performs a gravity model,
    and writes the output to `demand.csv`.

    :param node_file: Path to the node.csv file
    :param accessibility_matrix_file: Path to the accessibility_matrix.csv file
    :param output_demand_file: Path to the output demand.csv file
    :param beta: Friction factor for the gravity model
    :param model_type: Type of gravity model: "singly" or "doubly"
    """
    # Step 1: Read node.csv
    node_df = pd.read_csv(node_file)

    # Filter for zones with zone_id > 0 (activity zones)
    zones_df = node_df[node_df['zone_id'] > 0]

    # Create dictionaries for production and attraction values by zone_id
    production = zones_df.set_index('zone_id')['production'].to_dict()
    attraction = zones_df.set_index('zone_id')['attraction'].to_dict()

    # Step 2: Read accessibility_matrix.csv
    accessibility_df = pd.read_csv(accessibility_matrix_file)

    # Create empty list for storing demand data
    demand_data = []

    # Step 3: Prepare the model based on the type
    if model_type == "doubly":
        # Compute balancing factors for doubly constrained model
        B = {zone: 1.0 for zone in production.keys()}  # Initialize balancing factor B
        tolerance = 1e-6  # Tolerance for convergence
        max_iterations = 1000  # Maximum iterations
        iterations = 0
        converged = False

        # Iteratively compute balancing factors
        while not converged and iterations < max_iterations:
            converged = True
            A = {zone: 1.0 / sum(attraction[dest] * math.exp(-beta * accessibility_df[
                (accessibility_df["o_zone_id"] == zone) & (accessibility_df["d_zone_id"] == dest)
            ]["cost"].values[0]) * B[dest] for dest in attraction.keys() if dest != zone) for zone in production.keys()}

            for zone in production.keys():
                new_B = 1.0 / sum(production[orig] * A[orig] * math.exp(-beta * accessibility_df[
                    (accessibility_df["d_zone_id"] == orig) & (accessibility_df["d_zone_id"] == zone)
                ]["cost"].values[0]) for orig in production.keys() if orig != zone)
                if abs(new_B - B[zone]) > tolerance:
                    converged = False
                B[zone] = new_B

            iterations += 1

        if iterations == max_iterations:
            print("Warning: Balancing factors did not converge.")

    # Step 4: Calculate trip volumes using the chosen gravity model
    for _, row in accessibility_df.iterrows():
        from_zone = row['o_zone_id']
        to_zone = row['d_zone_id']
        cost = row['cost']

        # Skip if either zone has zero production/attraction
        if from_zone not in production or to_zone not in attraction:
            continue

        # Gravity model formula
        if model_type == "singly":
            Tij = production[from_zone] * attraction[to_zone] * math.exp(-beta * cost)
        elif model_type == "doubly":
            Tij = production[from_zone] * A[from_zone] * attraction[to_zone] * math.exp(-beta * cost) * B[to_zone]
        else:
            raise ValueError("Invalid model_type. Choose 'singly' or 'doubly'.")

        # Append to demand data
        demand_data.append({
            'o_zone_id': from_zone,
            'd_zone_id': to_zone,
            'volume': Tij
        })

    # Step 5: Write demand.csv
    demand_df = pd.DataFrame(demand_data)
    demand_df.to_csv(output_demand_file, index=False)

    print(f"Demand matrix has been successfully written to {output_demand_file}.")

           
class ModeType:
    """Represents a mode type with its associated attributes."""
    def __init__(self, mode_type, vot, pce, occ, dedicated_shortest_path, demand_file):
        self.mode_type = mode_type  # Mode type as a string
        self.vot = vot  # Value of time
        self.pce = pce  # Passenger car equivalent
        self.occ = occ  # Occupancy
        self.dedicated_shortest_path = dedicated_shortest_path  # 1 if dedicated, 0 otherwise
        self.demand_file = demand_file  # Path to the demand file


# Equivalent to the array g_mode_type_vector[MAX_MODE_TYPES]

g_mode_type_vector = [None] * MAX_MODE_TYPES  # Create a list for mode types

class LinkRecord:
    def __init__(self):
        self.internal_from_node_id = 0
        self.internal_to_node_id = 0
        self.link_id = 0
        self.link_type = 1
        self.external_from_node_id = 0
        self.external_to_node_id = 0

        self.Lane_Capacity = 0.0
        self.Link_Capacity = 0.0
        self.lanes = 1.0
        self.FreeTravelTime = 0.0
        self.free_speed = 10.0
        self.Cutoff_Speed = 0.0

        self.VDF_Alpha = 0.15
        self.VDF_Beta = 4
        self.VDF_plf = 1.0
        self.Q_cd = 1.0
        self.Q_n = 1.0
        self.Q_cp = 0.28125
        self.Q_s = 4.0

        self.length = 0.0
        self.Speed = 0.0
        self.allowed_uses = "all"

        self.mode_allowed_use = [0] * MAX_MODE_TYPES
        self.mode_MainVolume = [0.0] * MAX_MODE_TYPES

        self.mode_SubVolume = [0.0] * MAX_MODE_TYPES
        self.mode_SDVolume = [0.0] * MAX_MODE_TYPES

        self.mode_Toll = [0.0] * MAX_MODE_TYPES
        self.mode_AdditionalCost = [0.0] * MAX_MODE_TYPES

        self.travel_time = 0.0
        self.BPR_TT = 0.0
        self.QVDF_TT = 0.0

        self.GenCost = 0.0
        self.GenCostDer = 0.0
        self.Ref_volume = 0.0
        self.Base_volume = 0.0
        self.Obs_volume = -1.0
        self.geometry = ""

    def setup(self, num_of_modes):
        self.link_type = 1
        self.VDF_Alpha = 0.15
        self.VDF_Beta = 4
        self.VDF_plf = 1.0
        self.Q_cd = 1.0
        self.Q_n = 1.0
        self.Q_cp = 0.28125
        self.Q_s = 4.0
        self.travel_time = 0.0
        self.BPR_TT = 0.0
        self.QVDF_TT = 0.0
        self.Ref_volume = 0.0
        self.Base_volume = 0.0
        self.Obs_volume = -1.0

def alloc_1d(size, default_value=0.0):
    """
    Allocate a 1D array initialized with a default value.

    :param size: Size of the array.
    :param default_value: Value to initialize the array with.
    :return: NumPy 1D array.
    """
    return np.full(size + 1, default_value, dtype=float)  # +1 to mimic C++ behavior

def alloc_1d_int(size, default_value=0.0):
    """
    Allocate a 1D array initialized with a default value.

    :param size: Size of the array.
    :param default_value: Value to initialize the array with.
    :return: NumPy 1D array.
    """
    return np.full(size + 1, default_value, dtype=int)  # +1 to mimic C++ behavior


def alloc_2d(rows, cols, default_value=0.0):
    """
    Allocate a 2D array initialized with a default value.

    :param rows: Number of rows in the array.
    :param cols: Number of columns in the array.
    :param default_value: Value to initialize the array with.
    :return: NumPy 2D array.
    """
    return np.full((rows + 1, cols + 1), default_value, dtype=float)  # +1 to mimic C++ behavior
def alloc_2d_int(rows, cols, default_value=0.0):
    """
    Allocate a 2D array initialized with a default value.

    :param rows: Number of rows in the array.
    :param cols: Number of columns in the array.
    :param default_value: Value to initialize the array with.
    :return: NumPy 2D array.
    """
    return np.full((rows + 1, cols + 1), default_value, dtype=int)  # +1 to mimic C++ behavior
def alloc_3d(shape, default_value=0.0):
    """
    Allocate a 3D array with the given shape, initialized with a default value.

    :param shape: Tuple specifying the dimensions (e.g., (dim1, dim2, dim3)).
    :param default_value: Value to initialize the array with.
    :return: NumPy 3D array.
    """
    return np.full(shape, default_value, dtype=float)

def alloc_3d_int(shape, default_value=0.0):
    """
    Allocate a 3D array with the given shape, initialized with a default value.

    :param shape: Tuple specifying the dimensions (e.g., (dim1, dim2, dim3)).
    :param default_value: Value to initialize the array with.
    :return: NumPy 3D array.
    """
    return np.full(shape, default_value, dtype=int)


def initialize_g_link_indices():
    """Initialize the global 5D list for link sequences."""
    global g_link_indices, max_routes, number_of_zones, number_of_modes

    start_time = time.time()

    # Create a 5D list filled with empty lists
    g_link_indices = [
        [
            [
                [[] for _ in range(max_routes + 1)]
                for _ in range(number_of_zones + 1)
            ]
            for _ in range(number_of_zones + 1)
        ]
        for _ in range(number_of_modes + 1)
    ]

    end_time = time.time()
    duration = end_time - start_time

    hours, rem = divmod(duration, 3600)
    minutes, seconds = divmod(rem, 60)
    milliseconds = (seconds - int(seconds)) * 1000
    seconds = int(seconds)

    print(f"Memory creation time for 5D link path matrix: {int(hours)} hours {int(minutes)} minutes {seconds} seconds {int(milliseconds)} ms")

def add_link_sequence(m, orig, dest, route_id, link_ids):
    """Add a link sequence for a specific mode, origin, destination, and route."""
    global g_link_indices

    if not g_link_indices:
        return

    # Ensure we are within bounds before adding the link sequence
    if (
        0 <= m < len(g_link_indices) and
        0 <= orig < len(g_link_indices[m]) and
        0 <= dest < len(g_link_indices[m][orig]) and
        0 <= route_id < len(g_link_indices[m][orig][dest])
    ):
        g_link_indices[m][orig][dest][route_id] = link_ids
    else:
        print("Error: Invalid indices for adding link sequence.")

                            
def process_node_file(file_name, log_file=None):
    """
    Reads node.csv and processes its data.
    
    :param file_name: Path to the node.csv file.
    :param log_file: Path to the log file, if any.
    :return: A tuple with the number of nodes, number of zones, 
             first through node, and mappings.
    """
    # Initialize variables
    global g_map_external_node_id_2_node_seq_no, g_map_external_node_id_2_node_seq_no, number_of_zones, number_of_nodes, first_thru_node


    l_first_thru_node = first_thru_node
    tap_log_enabled = log_file is not None

    # Open the log file if provided
    logfile = open(log_file, 'w') if log_file else None

    # Open and read the CSV file
    with open(file_name, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        for row in csv_reader:
            # Read node_id and zone_id
            node_id = int(row.get("node_id", 0))
            # Handle missing or empty zone_id values
            zone_id_value = row.get("zone_id", "")
            zone_id = int(zone_id_value) if zone_id_value.isdigit() else 0

            # Error check: zone_id should match node_id
            if zone_id >= 1 and zone_id != node_id:
                print(f"Error: zone_id should be the same as node_id but zone_id = {zone_id}, node_id = {node_id}")

            # Update mappings
            g_map_node_seq_no_2_external_node_id[number_of_nodes + 1] = node_id
            g_map_external_node_id_2_node_seq_no[node_id] = number_of_nodes + 1  # Sequential number starts from 1


            g_map_external_node_id_2_zone_id[node_id] = zone_id

            # Update number_of_zones
            if zone_id >= 1 and zone_id > number_of_zones:
                number_of_zones = zone_id

            # Set first through node if not initialized
            if zone_id == 0 and l_first_thru_node == first_thru_node:  # Not initialized
                l_first_thru_node = number_of_nodes + 1

            # Log data if tap_log_file is enabled
            if tap_log_enabled:
                logfile.write(f"node_id = {node_id}, node_seq_no = {g_map_external_node_id_2_node_seq_no[node_id]}\n")

            # Increment node count
            number_of_nodes += 1

    # Close the log file if it was opened
    if logfile:
        logfile.close()

    print(g_map_node_seq_no_2_external_node_id)
            
    return



{}

def get_number_of_links_from_link_file(file_name):
    """
    Reads a CSV file and counts the number of links based on its records.
    
    :param file_name: Path to the link.csv file.
    :return: The total number of links in the file.
    """
    global number_of_links
    number_of_links = 0

    # Open and read the CSV file
    with open(file_name, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        for row in csv_reader:
            # Read link_id and from_node_id (values are not used in counting)
            link_id = int(row.get("link_id", 0))
            from_node_id = int(row.get("from_node_id", 0))
            
            # Increment the link count
            number_of_links += 1

    return number_of_links


def read_links(file_name, use_us_standard):
    global number_of_links, number_of_modes, g_map_external_node_id_2_node_seq_no, g_link_vector
    g_tap_log_file=0
    logfile_path=None
    #g_link_vector = defaultdict(LinkRecord)  # Dictionary to store links with their index as keys
    total_base_link_volume = 0
    k = 1  # Link index starts from 1

    with open(file_name, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)

        for row in csv_reader:
            link = LinkRecord()
            link.setup(number_of_modes)

            # Read basic properties
            link.external_from_node_id = int(row["from_node_id"])
            link.external_to_node_id = int(row["to_node_id"])
            link.link_id = int(row["link_id"])
            link.link_type = int(row["link_type"])

            # Map internal node IDs
            link.internal_from_node_id = g_map_external_node_id_2_node_seq_no.get(link.external_from_node_id, 0)
            if link.internal_from_node_id == 0:
                print(f"Error in from_node_id = {link.external_from_node_id} for link_id = {link.link_id}")
                continue

            link.internal_to_node_id = g_map_external_node_id_2_node_seq_no.get(link.external_to_node_id, 0)
            if link.internal_to_node_id == 0:
                print(f"Error in to_node_id = {link.external_to_node_id} for link_id = {link.link_id}")
                continue

            #  Set use_us_standard to True for US Standard (miles, mph), False for International (meters, kmph)


            
            link.Ref_volume = float(row.get("ref_volume", 0.0))
            link.length = float(row.get("length", 0.0))

            if "lanes" in row:
                link.lanes = int(row["lanes"])
            if "capacity" in row:
                link.Lane_Capacity = float(row["capacity"])
                link.Link_Capacity = link.lanes * link.Lane_Capacity

            link.free_speed = float(row.get("free_speed", 10.0))
            
            # Read additional properties
            # Determine the unit for length based on the standard
            if use_us_standard:
                # Assume length is in miles, speed in mph 
                pass  # No conversion needed
            else:
                # Length is in meters (International Standard), speed in kmph
                link.length /= 1609.34    # convert to mile
                link.free_speed  /=1.609  # convert to mile per hour

            
            link.FreeTravelTime = link.length / link.free_speed * 60.0

            debug = 1 
             # Debugging output
            if debug and k < 3: 
                print("Processing link {k} properties:")
                print(f"  Length: {link.length }")
                print(f"  Lanes: {link.lanes}")
                print(f"  Lane Capacity: {link.Lane_Capacity}")
                print(f"  Link Capacity: {link.Link_Capacity}")
                print(f"  Free Speed: {link.free_speed}")
                print(f"  Free Travel Time: {link.FreeTravelTime}")
        
            # # Handle mode-based properties
            # for m in range(1, number_of_modes + 1):
            #     mode_field = f"base_vol_{g_mode_type_vector[m]['mode_type']}"

            #     toll_field = f"toll_{g_mode_type_vector[m]['mode_type']}"
            #     if toll_field in row:
            #         link.mode_Toll[m] = float(row.get(toll_field, 0.0))
            #         link.mode_AdditionalCost[m] = link.mode_Toll[m] / g_mode_type_vector[m]["vot"] * 60.0

            # Final processing
            link.BoverC = link.VDF_Alpha / pow(link.Link_Capacity, link.VDF_Beta) if link.Link_Capacity > 0 else 0.0

            g_link_vector[k] = link
            k += 1

    print(f"total_base_link_volume = {total_base_link_volume}")
    baselinkvolume_loaded_flag = 1 if total_base_link_volume > 0 else 0

    return

def find_links_to():
    """
    Finds the sorted list of links directed to each node.

    :param no_nodes: Number of nodes.
    :param number_of_links: Number of links.
    :param links: List of Link objects indexed by 1-based index.
    :return: A dictionary where keys are node IDs, and values are sorted lists of link IDs.
    """
    global links_to, g_link_vector, number_of_links, number_of_nodes
    # Initialize a dictionary with empty lists for each node
    links_to = {node: [] for node in range(1, number_of_nodes + 1)}

    # Populate the links_to dictionary
    for k in range(1, number_of_links + 1):
        internal_to_node_id = g_link_vector[k].internal_to_node_id
        links_to[internal_to_node_id].append(k)

    # Sort the lists for each node
    for node in links_to:
        links_to[node].sort()

    return links_to

def init_link_pointers(links_file_name,  g_tap_log_file=0, logfile_path=None, use_us_standard = False):
    """
    Initializes pointers to the first and last link originating from each node.

    :param links_file_name: Name of the links file (for warnings and messages).
    :param no_nodes: Total number of nodes.
    :param number_of_links: Total number of links.
    :param links: List of Link objects indexed by 1-based index.
    :param g_map_node_seq_no_2_external_node_id: Map of internal to external node IDs.
    :param g_tap_log_file: Flag for logging output (1 = enabled, 0 = disabled).
    :param logfile_path: Path to the log file (if logging is enabled).
    :return: Tuple of FirstLinkFrom and LastLinkFrom lists.
    """
    # Initialize FirstLinkFrom and LastLinkFrom lists
    global FirstLinkFrom, LastLinkFrom , number_of_nodes, number_of_links
    FirstLinkFrom = [0] * (number_of_nodes + 1)  # 1-based indexing
    LastLinkFrom = [-1] * (number_of_nodes + 1)  # -1 indicates no links

    FirstLinkFrom[1] = 1
    Node = 1

    for k in range(1, number_of_links + 1):
        internal_from_node_id = g_link_vector[k].internal_from_node_id

        if internal_from_node_id == Node:
            continue
        elif internal_from_node_id >= Node + 1:
            LastLinkFrom[Node] = k - 1
            Node = internal_from_node_id
            FirstLinkFrom[Node] = k
        elif internal_from_node_id < Node:
            raise ValueError(f"Sort error in link file '{links_file_name}': "
                             f"a link from node {internal_from_node_id} was found after a link from node {Node}.")
        elif internal_from_node_id > Node + 1:
            # Handle nodes with no links
            LastLinkFrom[Node] = k - 1
            for Node in range(Node + 1, internal_from_node_id):
                FirstLinkFrom[Node] = 0
                LastLinkFrom[Node] = -1
            FirstLinkFrom[Node] = k

    if Node == number_of_nodes:
        LastLinkFrom[Node] = number_of_links
    else:
        # Handle remaining nodes with no links
        LastLinkFrom[Node] = number_of_links - 1
        for Node in range(Node + 1, no_nodes + 1):
            FirstLinkFrom[Node] = 0
            LastLinkFrom[Node] = -1

    # Optional logging
    if g_tap_log_file == 1 and logfile_path:
        with open(logfile_path, 'w') as logfile:
            for Node in range(1, number_of_nodes + 1):
                logfile.write(f"node_id = {g_map_node_seq_no_2_external_node_id.get(Node, Node)}, "
                              f"FirstLinkFrom = {FirstLinkFrom[Node]}, "
                              f"LastLinkFrom = {LastLinkFrom[Node]} \n")

    return




def init_links(links_file_name="link.csv", g_tap_log_file=0, logfile_path=None, use_us_standard= False):
    """
    Initializes link data by reading the links file, setting up the 'LinksTo' structure, and initializing link pointers.

    :param links_file_name: Name of the links file.

 
    :param g_tap_log_file: Flag for logging output (1 = enabled, 0 = disabled).
    :param logfile_path: Path to the log file (if logging is enabled).
    :return: Tuple containing links, LinksTo, FirstLinkFrom, and LastLinkFrom.
    """
    global g_link_vector, links_to, FirstLinkFrom, LastLinkFrom, number_of_links, number_of_modes, number_of_nodes, g_map_node_seq_no_2_external_node_id
    # Initialize the links structure
    g_link_vector = [None] * (number_of_links + 1)  # +1 for 1-based indexing

    # Step 1: Read links from file
    read_links(links_file_name, use_us_standard)

    print(f" after reading link.csv, Length of Link: {len(g_link_vector)}")


    # Step 2: Find 'LinksTo' structure
    links_to = find_links_to()

    # Step 3: Initialize link pointers
    init_link_pointers(links_file_name, g_tap_log_file, logfile_path,use_us_standard)

    return



def sum_od_table():
    """
    Calculate the sum of the OD table and populate total_o_table.
   
    """
    global md_od_flow, total_o_flow, number_of_zones, number_of_modes
    
    total_sum  = 0

    for m in range(1, number_of_modes+1):
        for orig in range(1, number_of_zones + 1):
            for dest in range(1, number_of_zones + 1):
                total_o_flow[orig] += md_od_flow[m, orig, dest]
                total_sum += md_od_flow[m, orig, dest]

    # Print the OD sum volume for each zone
    print("OD Sum Volume for Each Zone:")
    for orig in range(1, number_of_zones + 1):
        # Sum all destination flows for a given origin and mode
        zone_sum = md_od_flow[1, orig, :].sum()  # Sums all destinations for mode 1 and origin `orig`
        print(f"Zone {orig}: {zone_sum:.4f}")
            
    # Print the total OD volume
    print(f"Total OD Volume Across All Zones: {total_sum:.4f}")
    
    return total_sum

def read_od_flow():
    """
    Read the OD flow and initialize required data structures.
    """
    global number_of_modes, number_of_zones, g_link_vector,  md_od_flow, md_route_cost, total_o_flow, g_zone_outbound_link_size
    od_table_shape = (number_of_modes+1, number_of_zones + 1, number_of_zones + 1)
    md_od_flow = alloc_3d(od_table_shape)
    md_route_cost = alloc_3d(od_table_shape)
    total_o_flow = alloc_1d(number_of_zones + 1)

    m = 1 
    for orig in range(1, number_of_zones + 1):
        total_o_flow [orig] = 0
        for dest in range(1, number_of_zones + 1):
            md_od_flow [m][orig][dest] = 0 
            md_route_cost [m][orig][dest] = BIGM           


    # Populate MDODflow
    read_od_table()

    # Calculate the total OD flow
    real_total = sum_od_table()

    # Calculate outbound link size
    g_zone_outbound_link_size = alloc_1d(number_of_zones + 1, default_value=0)
    for k in range(1, number_of_links + 1):
        if  g_link_vector[k].external_from_node_id <= number_of_zones:
            g_zone_outbound_link_size[g_link_vector[k].external_from_node_id ] += 1

    # Error checking
    for z in range(1, number_of_zones + 1):
        if g_zone_outbound_link_size[z] == 0 and total_o_flow[z] > 0.01:
            raise ValueError(f"No outbound link from zone {z} with positive demand {total_o_flow[z]}")

    return

def read_od_table():
    """
    Read OD table and populate OD and DiffOD tables.
    """
    global number_of_zones,  md_od_flow
    
    m = 1 
    demand_file = 'demand.csv'
    print(f"Reading demand file: {demand_file}")

    line_count = 0

    try:
        with open(demand_file, 'r', encoding='utf-8-sig') as file:
            reader = csv.DictReader(file)
            for row in reader:
                o_zone_id = int(row["o_zone_id"])
                d_zone_id = int(row["d_zone_id"])
                volume = float(row["volume"])

                if o_zone_id > number_of_zones or d_zone_id > number_of_zones:
                    raise ValueError(f"Invalid zone ID in {demand_file}: {o_zone_id}, {d_zone_id}")

                if line_count <= 3:
                    print(f"o_zone_id: {o_zone_id}, d_zone_id: {d_zone_id}, volume: {volume:.4f}")
                    line_count = line_count + 1


                md_od_flow[m, o_zone_id, d_zone_id] = volume
                print(md_od_flow[m, o_zone_id, d_zone_id] )

    except FileNotFoundError:
        print(f"Error: File {demand_file} not found.")
        return
# Constants
INVALID = -1
BIGM = 9999999
WAS_IN_QUEUE = -7

def minpath(mode, orig, pred_link, cost_to):
    """
    Computes the shortest path from an origin node using a modified label-setting algorithm.
    """
    global FirstLinkFrom, LastLinkFrom, g_map_external_node_id_2_node_seq_no, number_of_zones, number_of_nodes, first_thru_node
    cost_to.fill(BIGM)
    pred_link.fill(INVALID)
    
    # Initialize variables
    queue_next = alloc_1d_int(number_of_nodes+1)
    
    now = g_map_external_node_id_2_node_seq_no[orig]
    internal_node_id_for_origin_zone = now
   
    pred_link[now] = INVALID
    cost_to[now] = 0.0

    scan_list = []
    scan_list.append(now)
    return2q_count = 0
    
    debug  = 0 

    
    if debug:
        print(f"Starting minpath computation for mode {mode}, origin {orig}")
        print(f"Initial node: {now}, queue initialized")

    while scan_list:
        now = scan_list.pop(0)  # Remove the first element from the scan list
        if now >= first_thru_node or now == internal_node_id_for_origin_zone:
            if debug:
                print(f"Processing node {now}...")
                
            for k in range(FirstLinkFrom[now], LastLinkFrom[now] + 1):
            #    if g_link_vector[k]["mode_allowed_use"][mode] == 0:
            #        continue

                new_node = g_link_vector[k].internal_to_node_id
                new_cost = cost_to[now] + g_link_vector[k].travel_time

                if debug:
                    print(f"Checking link {k}: new_node={new_node}, g_link_vector[k].length  = {g_link_vector[k].length}, new_cost={new_cost:.4f}, "
                          f"current_cost={cost_to[new_node]:.4f}")

                if cost_to[new_node] > new_cost:
                    
                    if debug:
                        print(f"Updated cost for node {new_node}: {new_cost:.4f}")
                        print(f"Predecessor for node {new_node}: link {k}")
    
                    cost_to[new_node] = new_cost
                    pred_link[new_node] = k

 # Add the node to the scan list if it's not already there
                    if new_node not in scan_list:
                        scan_list.append(new_node)
                        if debug:
                            print(f"    Node {new_node} added to scan list")

    if debug:
        print("Updated cost_to array:")
        for idx, cost in enumerate(cost_to):
            print(f"  Node {idx}: Cost = {cost:.4f}")
    
        print("Updated pred_link array:")
        for idx, pred in enumerate(pred_link):
            print(f"  Node {idx}: Predecessor Link = {pred}")
    
    if debug:
        print("Finished minpath computation")
        print(f"Total nodes returned to queue: {return2q_count}")
            
        return return2q_count

def find_min_cost_routes(min_path_pred_link):
    """
    Finds the minimum cost routes for all origins and modes.
    """
    global g_map_external_node_id_2_node_seq_no, g_map_external_node_id_2_zone_id, total_o_flow, md_route_cost, md_od_flow, g_zone_outbound_link_size, number_of_zones, number_of_nodes
    cost_to = np.full((number_of_zones + 1, number_of_nodes + 1), BIGM, dtype=float)

    system_least_travel_time  = 0 
    m = 1 
    
    
    debug = 0
    
    for orig in range(1, number_of_zones+1):

        if orig not in g_map_external_node_id_2_zone_id:
            continue
        
        if g_map_external_node_id_2_zone_id[orig] <1:
            continue
    

        
        minpath(m, orig, min_path_pred_link[m][orig], cost_to[orig])

        if md_route_cost is not None:
            for dest in range(1, number_of_zones + 1):
               if dest not in g_map_external_node_id_2_zone_id:
                    continue
        
               md_route_cost[m][orig][dest] = BIGM

               internal_node_id_for_destination_zone = g_map_external_node_id_2_node_seq_no[dest]
               if cost_to[orig][internal_node_id_for_destination_zone] <= BIGM - 1:
                    md_route_cost[m][orig][dest] = cost_to[orig][internal_node_id_for_destination_zone]
                    if md_od_flow[m][orig][dest] > 1e-6:
                        system_least_travel_time += (
                            md_route_cost[m][orig][dest] * md_od_flow[m][orig][dest]
                        )

                    if debug:
                        print(f"Processing mode {m}, origin {orig}")
                        print(f"    Destination {dest}:")
                        print(f"      Cost to Destination: {cost_to[orig][internal_node_id_for_destination_zone]:.4f}")
                        print(f"      Updated md_route_cost[{m}][{orig}][{dest}] = {md_route_cost[m][orig][dest]:.4f}")
                        print(f"      Updated System Least Travel Time: {system_least_travel_time:.4f}")
   
    return system_least_travel_time
 
def link_travel_time(k, volume):
    global g_link_vector, demand_period_starting_hours, demand_period_ending_hours
   
    print(f"Processing Link {k}")
    
    if k >= len(g_link_vector):
        print(f"Invalid index: {k} exceeds Link size {len(Link)}")
       
    


    incoming_demand = (
        volume[k] 
        / max(0.01, g_link_vector[k].lanes)
        / max(0.001, demand_period_ending_hours - demand_period_starting_hours)
        / max(0.0001, g_link_vector[k].VDF_plf)
    )
    print(f"Link {k}: Incoming Demand = {incoming_demand}")

    g_link_vector[k].travel_time = (
        g_link_vector[k].FreeTravelTime 
        * (1.0 + g_link_vector[k].VDF_Alpha * (incoming_demand / max(0.1, g_link_vector[k].Link_Capacity)) ** g_link_vector[k].VDF_Beta)
    )

    print(f"Link {k}: Calculated Travel Time = {g_link_vector[k].travel_time}")
    
    if g_link_vector[k].travel_time < 0:
        g_link_vector[k].travel_time = 0
        print(f"Link {k}: Travel Time adjusted to 0 (was negative)")

    g_link_vector[k].BPR_TT = g_link_vector[k].travel_time
    print(f"Link {k}: BPR_TT = {g_link_vector[k].BPR_TT}")

    return g_link_vector[k].travel_time

def link_travel_time_integral(k, volume):
    global g_link_vector,demand_period_starting_hours, demand_period_ending_hours
    incoming_demand = (
        volume[k]
        / max(0.001, demand_period_ending_hours - demand_period_starting_hours)
        / max(0.0001, g_link_vector[k].VDF_plf)
    )
    integral = 0.0
    if g_link_vector[k].VDF_Beta >= 0.0:
        integral += (
            incoming_demand
            + (
                volume[k]
                * g_link_vector[k].FreeTravelTime
                * (
                    1.0
                    + (g_link_vector[k].BoverC / (g_link_vector[k].VDF_Beta + 1))
                    * (incoming_demand ** (g_link_vector[k].VDF_Beta + 1))
                )
            )
        )
    return integral


def link_travel_time_derivative(k, volume):
    global g_link_vector 
    if g_link_vector[k].VDF_Beta == 0.0:
        return 0.0
    else:
        return (
            g_link_vector[k].FreeTravelTime
            * g_link_vector[k].BoverC
            * g_link_vector[k].VDF_Beta
            * (volume[k] ** (g_link_vector[k].VDF_Beta - 1))
        )


def additional_cost(k, m):
    global g_link_vector    
    add_cost = g_link_vector[k].mode_Toll[m] / g_mode_type_vector[m].vot * 60.0
    return add_cost


def link_gen_cost(k, volume):
    global g_link_vector 
    return g_link_vector[k].mode_AdditionalCost[1] + link_travel_time(k, volume)


def link_cost_integral(k, volume):
    global g_link_vector 
    return g_link_vector[k].mode_AdditionalCost[1] * volume[k] + link_travel_time_integral(k, volume)


def link_gen_cost_derivative(k, volume):
    return link_travel_time_derivative(k, volume)

def update_link_cost(main_volume):
    global g_link_vector , number_of_links
    system_wide_travel_time = 0.0

    print(f"Length of Link: {len(g_link_vector)}")

    for k in range(1, number_of_links + 1):  # Adjust range as needed
        # if k >= len(g_link_vector):
        #     print(f"Invalid index: {k} exceeds Link size {len(Link)}")




        #print(f"Link {k}: calling link_travel_time")
        g_link_vector[k].travel_time = link_travel_time(k, main_volume)
        #print(f"Link {k}: Travel Time = {g_link_vector[k].travel_time}")

        g_link_vector[k].GenCost = link_gen_cost(k, main_volume)
        #print(f"Link {k}: Generalized Cost = {g_link_vector[k].GenCost}")

        contribution_to_travel_time = main_volume[k] * g_link_vector[k].travel_time
        system_wide_travel_time += contribution_to_travel_time
        #print(f"Link {k}: Contribution to System-Wide Travel Time = {contribution_to_travel_time}")

    print(f"Total System-Wide Travel Time = {system_wide_travel_time}")
    return system_wide_travel_time

def output_assessibility_matrix():
    global number_of_zones, number_of_modes,  g_zone_outbound_link_size, g_map_external_node_id_2_node_seq_no, md_route_cost

    # Output file name
    output_file = "accessibility_matrix.csv"

    # Open the CSV file for writing
    with open(output_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the header row
        csv_writer.writerow(["o_zone_id", "d_zone_id", "mode_id", "cost"])

        # Loop through all origin zones
        for Orig in range(1, number_of_zones + 1):
            
            if Orig not in g_map_external_node_id_2_zone_id:
                continue

           
            # Skip zones with no outbound links
            if g_zone_outbound_link_size[Orig] == 0:
                continue

            # Iterate over all destinations
            for Dest in range(1, number_of_zones + 1):
               
                if Dest not in g_map_external_node_id_2_zone_id:
                   continue
                
                if Dest == Orig:
                    continue  # Skip self-loops

                for m in range(1, number_of_modes + 1):  # Iterate over modes
                    # Skip if cost is not feasible (BIGM represents an infeasible cost)
                    if md_route_cost[m][Orig][Dest] >= BIGM - 1:
                         continue

                    # Write the valid record to the CSV
                    csv_writer.writerow([Orig, Dest, m, md_route_cost[m][Orig][Dest]])

    print(f"Accessibility matrix has been successfully written to {output_file}.")

def all_or_nothing_assign(assignment_iteration_no, main_volume, MDMinPathPredLink):
    global g_link_vector, total_o_flow, number_of_zones,  number_of_modes, g_zone_outbound_link_size, md_od_flow, md_route_cost, g_map_external_node_id_2_node_seq_no  
    print(f"All or nothing assignment, assignment_iteration_no = {assignment_iteration_no}")
    
    # Initialize ProcessorVolume and ProcessorModeVolume
    Volume = np.zeros((number_of_links + 1))
    ModeVolume = np.zeros((number_of_links + 1, number_of_modes + 1))
    
    for Orig in range(1, number_of_zones + 1):
        # Skip zones with no positive flow or outbound links
        if total_o_flow[Orig] < 0.00001 or g_zone_outbound_link_size[Orig] == 0:
            continue
        
        for m in range(1, number_of_modes + 1):
            for Dest in range(1, number_of_zones + 1):
                if Dest == Orig:
                    continue
                
                RouteFlow = md_od_flow[m][Orig][Dest]
                if RouteFlow == 0:
                    continue
                
                if md_route_cost[m][Orig][Dest] >= BIGM - 1:
                    continue
                
                CurrentNode = g_map_external_node_id_2_node_seq_no[Dest]
                internal_node_for_origin_node = g_map_external_node_id_2_node_seq_no[Orig]
                currentLinkSequence = []
                
                while CurrentNode != internal_node_for_origin_node:
                    k = MDMinPathPredLink[1][Orig][CurrentNode]
                    
                    if k <= 0 or k > number_of_links:
                        print(f"A problem in All_or_Nothing_Assign() Invalid pred for node seq no {CurrentNode} Orig zone = {Orig}")
                        break
                    
                    Volume[k] += RouteFlow * 1 # g_mode_type_vector[m].pce
                    ModeVolume[k][m] += RouteFlow
                    CurrentNode = g_link_vector[k].internal_from_node_id
                    
                    if CurrentNode <= 0 or CurrentNode > number_of_nodes:
                        print(f"A problem in All_or_Nothing_Assign() Invalid node seq no {CurrentNode} Orig zone = {Orig}")
                        break
                    
                    currentLinkSequence.append(k)
                
                # Store link sequence for this OD pair
                add_link_sequence(m, Orig, Dest, assignment_iteration_no, currentLinkSequence)
    
    # Update volumes based on iteration
    if assignment_iteration_no == 0:
        for k in range(1, number_of_links + 1):
            main_volume[k] = 0 
            for m in range(1, number_of_modes + 1):
                g_link_vector[k].mode_MainVolume[m] = ModeVolume[k][m]
                main_volume[k] = main_volume[k]  + ModeVolume[k][m]
                print(f"AssignIterations=0: Link {k}, Mode {m}, mode_MainVolume updated to {g_link_vector[k].mode_MainVolume[m]}")
    else:
        for k in range(1, number_of_links + 1):
            main_volume[k] = 0 
            for m in range(1, number_of_modes + 1):
                main_volume[k] = main_volume[k]  + ModeVolume[k][m]
                g_link_vector[k].mode_SubVolume[m] = ModeVolume[k][m]
                print(f"AssignIterations>0: Link {k}, Mode {m}, mode_SubVolume updated to {g_link_vector[k].mode_SubVolume[m]}")


def write_link_performance():
    global demand_period_starting_hours, demand_period_ending_hours, g_link_vector
    """
    Write link performance data to a CSV file.

    Parameters:
        iteration_no: Current iteration number.
        main_volume: List of main volumes on links.
        links: List of link objects with attributes.
        mode_type_vector: List of mode type objects with attributes.
        demand_period_starting_hours: Start hour of the demand period.
        demand_period_ending_hours: End hour of the demand period.
    """
    # Open files for writing and appending
    link_performance_filename = "link_performance.csv"

    #print("g_link_vector initialization:")
    #for link in g_link_vector:
    #    print(link)
    
    try:
        with open(link_performance_filename, mode='w', newline='') as link_file:
            # Write headers
            writer = csv.writer(link_file)
            headers = [ "link_id", "from_node_id", "to_node_id", "volume", "ref_volume",
                "base_volume", "obs_volume", "capacity", "D", "doc", "fftt", "travel_time",
                "VDF_alpha", "VDF_beta", "VDF_plf", "speed", "VMT", "VHT", "PMT", "PHT",
                "VHT_QVDF", "PHT_QVDF", "geometry"
            ]
            # Add mode-specific headers
            #headers.extend([f"mod_vol_{mode.mode_type}" for mode in mode_type_vector])
            headers.extend([
                "P", "t0", "t2", "t3", "vt2", "mu", "Q_gamma", "free_speed", "cutoff_speed",
                "congestion_ref_speed", "avg_queue_speed", "avg_QVDF_period_speed",
                "avg_QVDF_period_travel_time", "Severe_Congestion_P"
            ])
            # Add speed intervals
            headers.extend(
                [f"spd_{t // 60:02}:{t % 60:02}" for t in range(demand_period_starting_hours * 60, demand_period_ending_hours * 60, 5)]
            )
            writer.writerow(headers)

        # Append data
        with open(link_performance_filename, mode='a', newline='') as link_file:
            writer = csv.writer(link_file)

            for k, link in enumerate(g_link_vector):
                # Placeholder values (replace with actual calculations)
                if link is None:
                    print(f"Error: Link at index {k} is None.")
                    continue
            
                if not hasattr(link, 'mode_MainVolume'):
                    print(f"Error: Link at index {k} is missing 'mode_MainVolume'.")
                    continue
                P = t0 = t2 = t3 = vt2 = mu = Severe_Congestion_P = 0
                Q_gamma = congestion_ref_speed = avg_queue_speed = avg_QVDF_period_speed = 0
                IncomingDemand = DOC = 0

                # Call your Link_QueueVDF function here
                model_speed = [0] * 300  # Placeholder for model speed array

                # Calculate VMT, VHT, PMT, PHT, VHT_QVDF, PHT_QVDF
                VMT = VHT = PMT = PHT = VHT_QVDF = PHT_QVDF = 0
                m  = 1 # for m, mode in enumerate(mode_type_vector, start=1):
                occ  =1 
                
                VMT += link.mode_MainVolume[m] * link.length
                VHT += link.mode_MainVolume[m] * link.travel_time / 60.0
                PMT += link.mode_MainVolume[m] * occ * link.length
                PHT += link.mode_MainVolume[m] * occ * link.travel_time / 60.0
                VHT_QVDF += link.mode_MainVolume[m] * link.QVDF_TT / 60.0
                PHT_QVDF += link.mode_MainVolume[m] * occ * link.QVDF_TT / 60.0

                # Create a row of data
                row = [
                    link.link_id, link.external_from_node_id, link.external_to_node_id,
                    link.mode_MainVolume[m], link.Ref_volume, link.Base_volume, link.Obs_volume, link.Link_Capacity,
                    IncomingDemand, DOC, link.FreeTravelTime, link.travel_time, link.VDF_Alpha, link.VDF_Beta,
                    link.VDF_plf, link.length / max(link.travel_time / 60.0, 0.001),
                    link.travel_time - link.FreeTravelTime, VMT, VHT, PMT, PHT, VHT_QVDF, PHT_QVDF,
                    link.geometry
                ]

                # Add mode-specific data m = 1
                row.extend([link.mode_MainVolume[m]] )

                # Add additional parameters
                row.extend([
                    P, t0, t2, t3, vt2, mu, Q_gamma, link.free_speed, link.Cutoff_Speed,
                    congestion_ref_speed, avg_queue_speed, avg_QVDF_period_speed, link.QVDF_TT, Severe_Congestion_P
                ])

                # Add speed data
                row.extend([model_speed[t // 5] for t in range(demand_period_starting_hours * 60, demand_period_ending_hours * 60, 5)])

                # Write the row to the file
                writer.writerow(row)
    except Exception as e:
        print(f"Error: {e}")


def output_route_details(filename='route_assignment.csv'):
    global g_link_indices, g_link_vector
    """
    Write route details to a CSV file.

    Args:
        filename (str): The name of the output CSV file.
        g_link_indices (list): Nested list containing route information.
        link_data (list): List of link objects with attributes such as length, travel_time, etc.
        mode_type_vector (list): List of mode types with `mode_type` attributes.
    """
    # Open the file for writing
    with open(filename, 'w', newline='') as output_file:
        writer = csv.writer(output_file)

        # Write the CSV header in lowercase
        writer.writerow([
            "mode", "route_seq_id", "o_zone_id", "d_zone_id", "unique_route_id", 
            "node_sequence", "link_ids", "total_distance", "total_free_flow_travel_time", 
            "total_travel_time", "route_key"
        ])

        if not g_link_indices:
            return  # Exit if no route data is available

        unique_route_id = 1
        # Loop through the modes
        for m in range(1, len(g_link_indices)):
            for orig in range(1, len(g_link_indices[m])):
                for dest in range(1, len(g_link_indices[m][orig])):
                    unique_routes = {}


                    for route_id, route in enumerate(g_link_indices[m][orig][dest]):
                        if route:  # Check if the route is non-empty
                            total_distance = 0.0
                            total_free_flow_travel_time = 0.0
                            total_travel_time = 0.0
                            node_ids_str = ""
                            link_ids_str = ""

                            node_sum = 0  # Sum of node IDs
                            link_sum = 0  # Sum of link IDs

                            # Collect node IDs, link indices, and compute total metrics
                            for i in range(len(route) - 1, -1, -1):
                                k = route[i]
                                #print(f"route[{i}] = {k}, type: {type(k)}")  # Debugging print statement
    
                                # Append the from_node_id for each link
                                from_node_id = g_link_vector[k].external_from_node_id
                                node_ids_str += f"{from_node_id};"
                                node_sum += from_node_id

                                # Append the link index
                                link_ids_str += f"{k};"
                                link_sum += k

                                # Sum up total distance and travel times
                                total_distance += g_link_vector[k].length
                                total_free_flow_travel_time += g_link_vector[k].FreeTravelTime
                                total_travel_time += g_link_vector[k].travel_time
                                # Add the to_node_id for the last link
                                if i == 0:
                                    to_node_id = g_link_vector[k].external_to_node_id
                                    node_ids_str += f"{to_node_id}"
                                    node_sum += to_node_id

                            # Create a unique key based on node sum and link sum
                            route_key = f"{node_sum}_{link_sum}"

                            # Check for uniqueness
                            if route_key not in unique_routes:
                                # Mark the route as unique
                                unique_routes[route_key] = True

                                # Remove trailing semicolon from link_ids_str
                                if link_ids_str.endswith(";"):
                                    link_ids_str = link_ids_str[:-1]

                                # Write the route data to the CSV file
                                writer.writerow([
                                    m, route_id, orig, dest, 
                                    unique_route_id, node_ids_str, link_ids_str, 
                                    total_distance, total_free_flow_travel_time, 
                                    total_travel_time, route_key
                                ])

                                unique_route_id = unique_route_id + 1

    print(f"Output written to {filename}")

def volume_difference(volume1, volume2, difference):
    global number_of_links, number_of_modes, g_link_vector
    for k in range(1, number_of_links + 1):
        difference[k] = volume1[k] - volume2[k]
        print(f"Link {k}: Volume1 = {volume1[k]}, Volume2 = {volume2[k]}, Difference = {difference[k]}")
        for m in range(1, number_of_modes + 1):
            g_link_vector[k].mode_MainVolume[m] = g_link_vector[k].mode_SubVolume[m] -  g_link_vector[k].mode_MainVolume[m] 
            
def update_volume(main_volume, sd_volume, lambda_):
    global number_of_links, number_of_modes, g_link_vector
    """
    Update the main volume using the search direction and step size lambda.

    Parameters:
        main_volume: Current flow volumes on each link (list or numpy array).
        sd_volume: Search direction volumes (difference between AON assignment and current flows).
        lambda_: Step size for updating volumes.
        link_data: A dictionary or data structure holding link information, including mode_MainVolume and mode_SDVolume.
        number_of_links: Total number of links.
        number_of_modes: Total number of modes.
    """
    # Update MainVolume using Lambda * SDVolume
    for k in range(1, number_of_links + 1):
        main_volume[k] += lambda_ * sd_volume[k]
        m  =1 
        g_link_vector[k].mode_MainVolume[m] = main_volume[k] 
        print(f"SDVolume = {sd_volume[k]}, Lambda = {lambda_}, Updated MainVolume = {main_volume[k]}")

# =============================================================================
#     # Update mode_MainVolume using Lambda * mode_SDVolume for each link and mode
#     for k in range(1, number_of_links + 1):
#         for m in range(1, number_of_modes + 1):
#             g_link_vector[k].mode_MainVolume[m] += lambda_ * g_link_vector[k].mode_SDVolume [m]
# =============================================================================
        

def of_links_directional_derivative(main_volume, sd_volume, Lambda):
    global number_of_links
    OFscale = 1 
    volume = np.zeros(number_of_links + 1)
    link_cost_sum = 0

    for k in range(1, number_of_links + 1):
        volume[k] = main_volume[k] + Lambda * sd_volume[k]

    for k in range(1, number_of_links + 1):
        link_cost_sum += link_gen_cost(k, volume) * sd_volume[k]

    return link_cost_sum / OFscale

def links_sd_line_search(main_volume, sd_volume):
    min_iterations  = 5 
    max_iterations = 5 
    
    """
    Perform a line search using the bisection method to find the optimal step size lambda.
    
    Parameters:
        main_volume: Current flow volumes on each link (list or numpy array).
        sd_volume: Search direction volumes (difference between AON assignment and current flows).
        of_links_directional_derivative: Function to compute the directional derivative.
        min_iterations: Minimum iterations for the bisection method.
        max_iterations: Maximum iterations for additional convergence steps.
    
    Returns:
        Optimal step size (lambda).
    """
    lambdaleft = 0
    lambdaright = 1
    lambda_ = 0.5

    # Initial check at lambda = 0
    grad = of_links_directional_derivative(main_volume, sd_volume, 0.0)

    if grad >= 0:
        global LastLambda
        LastLambda = 0.0
        return 0.0

    # Check at lambda = 1
    grad = of_links_directional_derivative(main_volume, sd_volume, 1.0)
    if grad <= 0:
        LastLambda = 1.0
        return 1.0

    # Bisection method for line search within [0, 1]
    for n in range(1, min_iterations + 1):
        grad = of_links_directional_derivative(main_volume, sd_volume, lambda_)

        if grad <= 0.0:
            lambdaleft = lambda_
        else:
            lambdaright = lambda_

        lambda_ = 0.5 * (lambdaleft + lambdaright)

    # Additional iterations to ensure convergence, if necessary max_iterations is MAX_NO_BISECTITERATION
    while lambdaleft == 0 and n <= max_iterations:
        grad = of_links_directional_derivative(main_volume, sd_volume, lambda_)

        if grad <= 0.0:
            lambdaleft = lambda_
        else:
            lambdaright = lambda_

        lambda_ = 0.5 * (lambdaleft + lambdaright)
        n += 1

    global ActualIterations
    ActualIterations = n - 1
    LastLambda = lambdaleft
    return lambdaleft
 





def traffic_assignment(assign_iterations,use_us_standard):

    # Example Usage
    file_name = "node.csv"  # Replace with your CSV file path
    log_file = "logfile.txt"  # Optional log file path
    
    # Process the node file
    process_node_file(file_name, log_file=log_file)
    
    
    
    # Example Usage
    file_name = "link.csv"  # Replace with your CSV file path
    
    # Get the number of links
    number_of_links = get_number_of_links_from_link_file(file_name)
    
    # Output the result
    print(f"Number of Links: {number_of_links}")
    
    
    # Output results
    print(f"Number of Nodes: {number_of_nodes}")
    print(f"Number of Zones: {number_of_zones}")
    print(f"First Through Node: {first_thru_node}")
    #print(f"Node Sequence Map: {g_map_node_seq_no_2_external_node_id}")
    #print(f"External Node Map: {external_node_map}")
    
    
    # Open the summary log file for writing
    summary_log_file = open("summary_log_file.txt", "w")
    
    # Initialize variables
    MainVolume = None  # Placeholder for a double array
    SubVolume = None  # Placeholder for a double array
    SDVolume = None  # Placeholder for a double array
    Lambda = 0.0
    MDMinPathPredLink = None  # Placeholder for a 3D array
    
    init_links(
        links_file_name="link.csv",
        g_tap_log_file=1,
        logfile_path="logfile.txt",
        use_us_standard = use_us_standard
    )
    
    print(f" before link index, Length of Link: {len(g_link_vector)}")
        
    initialize_g_link_indices()
    
    print(f" after link index, Length of Link: {len(g_link_vector)}")
        
    # Print results
    print("Links:", g_link_vector)
    print("LinksTo:", links_to)
    print("FirstLinkFrom:", FirstLinkFrom)
    print("LastLinkFrom:", LastLinkFrom)
    
    
    read_od_flow() 
    
    # Allocate MDMinPathPredLink
    
    iteration_no = 0
    # Initialize arrays
    MainVolume = np.zeros(number_of_links + 1)
    SDVolume = np.zeros(number_of_links + 1)
    SubVolume = np.zeros(number_of_links + 1)
    MDMinPathPredLink = np.zeros((number_of_modes + 1, number_of_zones + 1, number_of_nodes + 1), dtype=int)
    
    # Record the start time
    start = datetime.now()
    
    # Assign base volume to main volume
    for k in range(1, number_of_links + 1):
        MainVolume[k] = 0
    
    # Set up the cost using FFTT
    system_wide_travel_time = update_link_cost(MainVolume)
    
    # # Define the file path
    # file_path = "link_performance.csv"
    
    # # Open the file in write mode
    # with open(file_path, 'w') as link_performance_file:
    #     # Write headers to the file
    #     link_performance_file.write(
    #         "iteration_no,link_id,from_node_id,to_node_id,volume,ref_volume,base_volume,obs_volume,"
    #         "capacity,D,doc,fftt,travel_time,VDF_alpha,VDF_beta,VDF_plf,speed,VMT,VHT,PMT,PHT,"
    #         "VHT_QVDF,PHT_QVDF,geometry\n"
    #     )
    
    
    
    # for m in range(1, number_of_modes + 1):
    #     link_performance_file.write(f"mod_vol_{g_mode_type_vector[m].mode_type},")
    # link_performance_file.write(
    #     "P,t0,t2,t3,vt2,mu,Q_gamma,free_speed,cutoff_speed,congestion_ref_speed,"
    #     "avg_queue_speed,avg_QVDF_period_speed,avg_QVDF_period_travel_time,"
    #     "Severe_Congestion_P,"
    #)
    
    # # Write time interval speeds
    # for t in range(demand_period_starting_hours * 60, demand_period_ending_hours * 60, 5):
    #     hour = t // 60
    #     minute = t % 60
    #     link_performance_file.write(f"spd_{hour:02d}:{minute:02d},")
    # link_performance_file.write("\n")
    # Calculate the system least travel time
    
    
    
    system_least_travel_time = find_min_cost_routes(MDMinPathPredLink)
    output_assessibility_matrix()
    all_or_nothing_assign(0, MainVolume,MDMinPathPredLink)
    
    system_wide_travel_time = update_link_cost(MainVolume)
    
    
    for iteration_no in range(1, assign_iterations):
        system_least_travel_time = find_min_cost_routes(MDMinPathPredLink)
        
           
        all_or_nothing_assign(iteration_no, SubVolume,MDMinPathPredLink)  # assign to the subvolume 
    
        volume_difference(SubVolume, MainVolume, SDVolume)
    
        lambda_ = links_sd_line_search(MainVolume, SDVolume)
    
        update_volume(MainVolume, SDVolume, lambda_)
    
        system_wide_travel_time = update_link_cost(MainVolume)
    
        gap = (system_wide_travel_time - system_least_travel_time) / max(0.1, system_least_travel_time) * 100
    
    #        print(f"Iter No = {iteration_no}, Lambda = {lambda_:.6f}, System VMT = {system_wide_travel_time:.1f}, Least TT = {system_least_travel_time:.1f}, Gap = {gap:.2f}%")
            
    output_route_details('route_assignment.csv')
    write_link_performance()
    # Output
    print(f"System Least Travel Time: {system_least_travel_time:.2f}")
 
# main program below   
#osm2gmns_network()
#trip_generation()   
use_us_standard = False
#traffic_assignment(0,use_us_standard) 
#perform_simple_logit_model_for_trip_distribution(0.1)

#use_us_standard = False
traffic_assignment(2,use_us_standard)