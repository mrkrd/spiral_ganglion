#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = "spiral_ganglion",
    version = "0.2",
    packages = find_packages(),
    package_data = {
        "spiral_ganglion": ["*.mod"]
    },

    author = "Marek Rudnicki",
    author_email = "marek.rudnicki@tum.de",
    description = "Spiral ganglion model.",
    license = "GPL",
)
