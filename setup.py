#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.rst') as file:
    long_description = file.read()

setup(
    name = "spiral_ganglion",
    version = "1",
    author = "Marek Rudnicki",
    author_email = "marek.rudnicki@tum.de",

    description = "Spiral ganglion model.",
    license = "GPLv3+",

    packages = find_packages(),
    package_data = {
        "spiral_ganglion": ["*.mod"]
    },
    long_description = long_description,
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
    ],

    platforms = ["Linux", "FreeBSD"],
    install_requires=["numpy", "neuron", "cochlea"],
)
