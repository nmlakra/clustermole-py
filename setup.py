from setuptools import find_packages, setup

setup(
    name="clustermolepy",  # Package name
    version="0.3.0",  # Version
    python_requries=">=python3.10",
    packages=find_packages(),  # Automatically finds all packages
    install_requires=[
        "pandas",
        "requests",
        "biomart",
        "anndata",
        "numpy",
    ],  # Dependencies
)
