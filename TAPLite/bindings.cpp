#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "TAPLite.h"

namespace py = pybind11;

PYBIND11_MODULE(TAPLite, m) {
    m.doc() = "Python bindings for TAPLite traffic assignment model";

    m.def("read_settings_file", &read_settings_file, "Reads the settings file for TAPLite");
    m.def("read_mode_type_file", &read_mode_type_file, "Reads the mode type file for TAPLite");
    m.def("get_number_of_nodes", &get_number_of_nodes_from_node_file,
          "Returns the number of nodes from the node file",
          py::arg("number_of_zones"), py::arg("first_thru_node"));
    m.def("get_number_of_links", &get_number_of_links_from_link_file,
          "Returns the number of links from the link file");

    m.def("initialize", &Init, "Initializes TAPLite",
          py::arg("number_of_modes"), py::arg("no_zones"));
    m.def("initialize_link_indices", &InitializeLinkIndices,
          "Initializes link indices",
          py::arg("number_of_modes"), py::arg("no_zones"), py::arg("total_assign_iterations"));

    m.def("find_min_cost_routes", &FindMinCostRoutes,
          "Finds minimum cost routes");
    m.def("all_or_nothing_assign", &All_or_Nothing_Assign,
          "Performs all-or-nothing assignment");
    m.def("update_link_cost", &UpdateLinkCost,
          "Updates the link costs");
    m.def("volume_difference", &VolumeDifference,
          "Computes volume differences");
    m.def("links_sd_line_search", &LinksSDLineSearch,
          "Performs line search optimization");
    m.def("update_volume", &UpdateVolume,
          "Updates volume values after line search");

    // Logging
    m.def("exit_message", &ExitMessage, "Logs an exit message");

    // Constants
    m.attr("BUFFERSIZE") = BUFFERSIZE;
    m.attr("MAX_NO_BISECT_ITERATION") = MAX_NO_BISECT_ITERATION;
}
