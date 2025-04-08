""" Find shortest path given a from node and a to node

Two path engines are provided:
1. C++ engine which is a special implementation of the deque implementation in
   C++ and built into path_engine.dll.
2. Python engine which provides three implementations: FIFO, Deque, and
   heap-Dijkstra. The default is deque.
"""


import collections
import ctypes
import heapq
import platform
import os
from time import time

__all__ = [
    'assignment',
    'simulation'
]

# current_os = platform.system()
# if current_os == "Darwin":
#     library_name = "DTALite_arm.dylib"
# elif current_os == "Windows":
#     library_name = "DTALite.dll"
# elif current_os == "Linux":
#     library_name = "DTALite.so"
# else:
#     raise OSError("Unsupported operating system")


current_os = platform.system()
arch = platform.machine()
if current_os == "Darwin":
    if arch == "arm64":
        library_name = "DTALite_arm.dylib"
    elif arch == "x86_64":
        library_name = "DTALite_x86.dylib"
    else:
        raise OSError(f"Unsupported macOS architecture: {arch}")
elif current_os == "Windows":
    library_name = "DTALite.dll"
elif current_os == "Linux":
    library_name = "DTALite.so"
else:
    raise OSError("Unsupported operating system")


dtalib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "bin", library_name))


def assignment():
    dtalib.DTA_AssignmentAPI()

def simulation():
    dtalib.DTA_SimulationAPI()

