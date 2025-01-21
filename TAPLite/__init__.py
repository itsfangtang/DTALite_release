# __init__.py

# Import the Python bindings exposed from bindings.cpp
from .TAPLite import (
    read_settings_file,
    read_mode_type_file,
    get_number_of_nodes,
    get_number_of_links,
    initialize,
    initialize_link_indices,
    find_min_cost_routes,
    all_or_nothing_assign,
    update_link_cost,
    volume_difference,
    links_sd_line_search,
    update_volume,
)

# Metadata for the package
__version__ = "0.1.0"  # Define the version of your package
__author__ = "Xuesong ZHou, Han Zheng"  # Replace with your name
__license__ = "MIT"  # Or the license you're using

# Example: Exported functions for quick access
__all__ = [
    "read_settings_file",
    "read_mode_type_file",
    "get_number_of_nodes",
    "get_number_of_links",
    "initialize",
    "initialize_link_indices",
    "find_min_cost_routes",
    "all_or_nothing_assign",
    "update_link_cost",
    "volume_difference",
    "links_sd_line_search",
    "update_volume",
]
