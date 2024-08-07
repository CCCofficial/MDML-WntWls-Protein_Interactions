from setuptools import setup, find_packages

setup(
    name='IBM',
    version='0.1.0',
    packages=find_packages(include=['IBM', 'IBM.*']),
    install_requires=['numpy>=1.14.5','pytest',"MDAnalysis","pandas"]
)
