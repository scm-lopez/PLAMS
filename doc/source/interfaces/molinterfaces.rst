Molecule interfaces
-------------------------

In this section we present modules that allow translating between PLAMS molecules and similar objects from other Python packages.
Thanks to that, various methods of manipulating your system's geometry from these packages can be used within PLAMS scripts.



.. _ASEInterface:

ASE
~~~~~~~~~~~

"The Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase>`_) is a set of tools and Python modules for setting up, manipulating, running, visualizing and analyzing atomistic simulations.
The PLAMS interface to ASE is limited to handling |Molecule| objects.
It features access to the ``ase.io`` module for reading/writing |Molecule| objects and two functions that translate PLAMS |Molecule| objects into ASE ``Atoms`` objects and vice versa.

Please refer to the ASE documentation to see what can be done with ASE ``Atoms`` and its I/O module.

.. automodule:: scm.plams.interfaces.molecule.ase



.. _RDKitmol:

RDKit
~~~~~~~~~~~

`RDKit <https://www.rdkit.org>`_ is a collection of cheminformatics and machine-learning software written in C++ and Python.
PLAMS interface to RDKit originates from `QMFlows <https://github.com/SCM-NV/qmflows>`_ project and features functions for generating PLAMS molecules from string representations (SMARTS, SMILES) as well as a handful of tools for dealing with proteins and PDB files.

.. automodule:: scm.plams.interfaces.molecule.rdkit

