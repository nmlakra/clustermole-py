from setuptools import setup, find_packages

setup(
    name="clustermole_py",          # Package name
    version="0.1.0",           # Version
    packages=find_packages(),  # Automatically finds all packages
    install_requires=[         # Dependencies
        "pandas", 
        "requests"
    ],
)

