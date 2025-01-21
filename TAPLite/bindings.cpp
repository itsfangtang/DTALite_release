#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For handling STL containers if needed
#include "TAPLite.h"       // Include your TAPLite header file

namespace py = pybind11;

PYBIND11_MODULE(TAPLite, m) {
    m.doc() = "Python bindings for TAPLite traffic assignment model";

    // Initialization and setup functions
    m.def("read_settings_file", &read_settings_file, "Reads the settings file for TAPLite");
    m.def("read_mode_type_file", &read_mode_type_file, "Reads the mode type file for TAPLite");
    m.def("get_number_of_nodes", &get_number_of_nodes_from_node_file,
          "Returns the number of nodes from the node file",
          py::arg("no_zones"), py::arg("FirstThruNode"));
    m.def("get_number_of_links", &get_number_of_links_from_link_file,
          "Returns the number of links from the link file");

    // Initialization
    m.def("initialize", &Init, "Initializes TAPLite with the number of modes and zones",
          py::arg("number_of_modes"), py::arg("no_zones"));
    m.def("initialize_link_indices", &InitializeLinkIndices,
          "Initializes link indices for TAPLite",
          py::arg("number_of_modes"), py::arg("no_zones"), py::arg("TotalAssignIterations"));

    // Core assignment functions
    m.def("find_min_cost_routes", &FindMinCostRoutes,
          "Finds minimum cost routes",
          py::arg("MDMinPathPredLink"));
    m.def("all_or_nothing_assign", &All_or_Nothing_Assign,
          "Performs all-or-nothing assignment",
          py::arg("iteration_no"), py::arg("MDDiffODflow"),
          py::arg("MDMinPathPredLink"), py::arg("MainVolume"));

    // Volume and cost functions
    m.def("update_link_cost", &UpdateLinkCost,
          "Updates link cost given a volume array",
          py::arg("MainVolume"));
    m.def("volume_difference", &VolumeDifference,
          "Computes the volume difference between subvolume and main volume",
          py::arg("SubVolume"), py::arg("MainVolume"), py::arg("SDVolume"));
    m.def("links_sd_line_search", &LinksSDLineSearch,
          "Performs line search to optimize link volumes",
          py::arg("MainVolume"), py::arg("SDVolume"));
    m.def("update_volume", &UpdateVolume,
          "Updates volume values after line search",
          py::arg("MainVolume"), py::arg("SDVolume"), py::arg("Lambda"));

    // Logging and utilities
    m.def("log_summary", [](const std::string &message) {
        if (summary_log_file) {
            fprintf(summary_log_file, "%s\n", message.c_str());
        }
    }, "Logs a summary message to the summary log file");

    m.def("log_iteration", [](int iteration_no, double system_tt, double least_tt, double gap) {
        if (summary_log_file) {
            fprintf(summary_log_file,
                    "Iteration %d: System TT = %.2f, Least TT = %.2f, Gap = %.2f%%\n",
                    iteration_no, system_tt, least_tt, gap);
        }
    }, "Logs an iteration summary to the log file");

    // Expose constants (ensure these are defined in TAPLite.h)
    m.attr("BUFFERSIZE") = BUFFERSIZE;
    m.attr("MAX_NO_BISECT_ITERATION") = MAX_NO_BISECT_ITERATION;
}
