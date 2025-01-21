import TAPLite

def main():
    # Open summary log file
    summary_log_file = open("summary_log_file.txt", "w")

    # Step 1: Read Input Files
    TAPLite.read_settings_file()
    TAPLite.read_mode_type_file()
    print("Settings and mode type files read successfully.")

    # Step 2: Initialize Network and Links
    no_zones = 100  # Example value, set according to your data
    first_thru_node = 1
    total_assign_iterations = 10
    number_of_modes = 3

    num_nodes = TAPLite.get_number_of_nodes(no_zones, first_thru_node)
    num_links = TAPLite.get_number_of_links()

    print(f"Number of Nodes: {num_nodes}, Number of Links: {num_links}")
    summary_log_file.write(f"Number of Nodes: {num_nodes}, Number of Links: {num_links}\n")

    TAPLite.initialize(number_of_modes, no_zones)
    TAPLite.initialize_link_indices(number_of_modes, no_zones, total_assign_iterations)

    # Step 3: Initialize Variables
    MainVolume = [0] * num_links
    SubVolume = [0] * num_links
    SDVolume = [0] * num_links
    MDMinPathPredLink = [[[0] * num_nodes for _ in range(no_zones)] for _ in range(number_of_modes)]

    # Step 4: Iterative Traffic Assignment
    for iteration_no in range(1, total_assign_iterations + 1):
        print(f"Starting iteration {iteration_no}...")

        # Find minimum cost routes
        least_travel_time = TAPLite.find_min_cost_routes(MDMinPathPredLink)

        # Perform all-or-nothing assignment
        TAPLite.all_or_nothing_assign(iteration_no, [0] * num_links, MDMinPathPredLink, MainVolume)

        # Update link costs
        system_wide_travel_time = TAPLite.update_link_cost(MainVolume)

        # Compute the gap
        gap = (system_wide_travel_time - least_travel_time) / max(0.1, least_travel_time) * 100
        print(f"Iteration {iteration_no}: Gap = {gap:.2f}%")
        summary_log_file.write(f"Iteration {iteration_no}: Gap = {gap:.2f}%\n")

    # Close summary log file
    summary_log_file.close()
    print("Traffic assignment completed.")

if __name__ == "__main__":
    main()
