from setuptools import setup, find_packages

setup(
    name="clustermolepy",  # Package name
    version="0.2.0",  # Version
    packages=find_packages(),  # Automatically finds all packages
    install_requires=["pandas", "requests", "biomart", "anndata", "numpy"],  # Dependencies
)
