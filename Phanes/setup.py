from setuptools import setup

setup(
    name="phanes",
    version="0.1",
    description='Biological Sequences Handling',
    url='https://github.com/Educorreia932/FCUP-BINF/Phanes',
    author="Bruno Vaz, Eduardo Correia, Filipe Justiça",
    packages=["phanes", "phanes.alignment", "phanes.sequence"],
    install_requires=[
        "ordered_set"
    ],
    keywords=["python", "bioninformatics", "dna", "protein", "sequence-alignment"],
)
