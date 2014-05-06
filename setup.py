#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = "spiral_ganglion",
    version = "0.1",
    packages = find_packages(),
    package_data = {
        "spiral_ganglion": ["*.mod"]
    },

    author = "Marek Rudnicki",
    author_email = "marek.rudnicki@tum.de",
    description = "Spiral Ganglion model in Python.",
    license = "GPL",
)
