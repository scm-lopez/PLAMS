Molecule
-------------------------

In this chapter the PLAMS module responsible for handling molecular geometries is presented.
Information about atomic coordinates can be read from (or written to) files of various types: ``xyz``, ``pdb``, ``mol`` or ``mol2``.
PLAMS not only extracts the relevant data from those files, but also tries to "understand" the structure of the underlying molecule in terms of atoms and bonds between them, allowing performing a variety of simple operations like moving or rotating some parts of the molecule, splitting it into multiple parts, merging two molecules etc.
Further file types can be read using the Atomistic Simulation Environment (ASE), see :ref:`ASEInterface`.
Classes defined in this module are |Molecule|, |Atom| and |Bond|.

.. note::

    For importing coordinates in other formats and performing more advanced operations you can use |Molecule| interfaces to :ref:`ASEInterface` and :ref:`RDKitmol`.



Molecule
~~~~~~~~~~~~~~~~~~~~~~~~~


.. autoclass :: scm.plams.mol.molecule.Molecule
    :exclude-members: __weakref__, __copy__


Atom labeling
+++++++++++++

This subsection describes API of ``identify`` module, which is used to assign unique names to atoms in a molecule. Unique atom names are used in |Molecule| labeling (:meth:`~scm.plams.mol.molecule.Molecule.label`) and in the method restoring the order of atoms (:meth:`~scm.plams.mol.molecule.Molecule.reorder`). All the functions below, apart from :func:`~scm.plams.mol.identify.label_atoms`, are for internal use and they are not visible in the main PLAMS namespace.


.. currentmodule:: scm.plams.mol.identify

.. autofunction:: label_atoms
.. autofunction:: find_permutation
.. autofunction:: molecule_name
.. autofunction:: initialize
.. autofunction:: iterate
.. autofunction:: clear
.. autofunction:: new_name
.. autofunction:: knock
.. autofunction:: twist
.. autofunction:: bend
.. autofunction:: unique_atoms


Atom
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: scm.plams.mol.atom.Atom
    :exclude-members: __weakref__


Bond
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass :: scm.plams.mol.bond.Bond
    :exclude-members: __weakref__

