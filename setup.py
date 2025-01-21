from skbuild import setup
from setuptools import find_packages
import os
import sys

setup(
    name="TAPLite",
    version="0.1.0",
    description="A Python binding for TAPLite, a traffic assignment problem solver",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="Dr. Xuesong (Simon) Zhou, Dr. Han Zheng",
    author_email="your_email@example.com",  # Replace with your email
    url="https://github.com/your_username/TAPLite",  # Replace with your repo URL
    license="MIT",
    packages=find_packages(),
    cmake_install_dir="TAPLite",
    cmake_args=[
        "-DCMAKE_BUILD_TYPE=Release",
        "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable)
    ],
    python_requires=">=3.7",
    install_requires=[
        "pybind11>=2.6.0",
        "scikit-build",
        "setuptools",
        "wheel",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
)
