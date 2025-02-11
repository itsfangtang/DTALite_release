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

current_os = platform.system()
if current_os == "Darwin":
    library_name = "DTALite_arm.dylib"
elif current_os == "Windows":
    library_name = "DTALite.dll"
elif current_os == "Linux":
    library_name = "DTALite.so"
else:
    raise OSError("Unsupported operating system")

dtalib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), library_name))


def assignment():
    dtalib.DTA_AssignmentAPI()

def simulation():
    dtalib.DTA_SimulationAPI()

