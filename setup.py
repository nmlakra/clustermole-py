from setuptools import find_packages, setup

setup(
    name="clustermolepy",  # Package name
    version="0.2.0",  # Version
    python_requies=">=python3.12",
    packages=find_packages(),  # Automatically finds all packages
    install_requires=[
        "pandas",
        "requests",
        "biomart",
        "anndata",
        "numpy",
    ],  # Dependencies
)
