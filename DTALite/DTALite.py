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
from os import path
from time import time

__all__ = [
    'assignment',
    'simulation'
]


_os = platform.system()
if _os.startswith('Windows'):
    _dll_file = path.join(path.dirname(__file__), 'bin/DTALite.dll')
elif _os.startswith('Linux'):
    _dll_file = path.join(path.dirname(__file__), 'bin/DTALite.so')
elif _os.startswith('Darwin'):
    # check CPU is Intel or Apple Silicon
    if platform.machine().startswith('x86_64'):
        _dll_file = path.join(path.dirname(__file__), 'bin/DTALite_x86.dylib')
    else:
        _dll_file = path.join(path.dirname(__file__), 'bin/DTALite_arm.dylib')
else:
    raise Exception('Please build the shared library compatible to your OS\
                    using source files in DTALite/core!')

_cdll = ctypes.cdll.LoadLibrary(_dll_file)


def assignment():
    _cdll.DTA_AssignmentAPI()

def simulation():
    _cdll.DTA_SimulationAPI()

