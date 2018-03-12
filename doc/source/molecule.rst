Molecule
-------------------------

.. currentmodule:: scm.plams.core.basemol

In this chapter the PLAMS module responsible for handling molecular geometries is presented. Information about atomic coordinates can be read from (or written to) files of various types: ``xyz``, ``pdb``, ``mol`` or ``mol2``. PLAMS not only extracts relevant data from those files, but also tries to "understand" the structure of the underlying molecule in terms of atoms and bonds between them, allowing you to perform a variety of simple operations like, for example, moving or rotating some parts of the molecule, splitting it into multiple parts, merging two molecules etc. Further file types can be read using the Atomistic Simulation Environment (ASE), see :meth:`read`.  Classes defined in this module are |Molecule|, |Atom| and |Bond|. They interact with each other to provide a basic set of functionalities for geometry handling.



Molecule
~~~~~~~~~~~~~~~~~~~~~~~~~


.. autoclass :: Molecule
    :exclude-members: __weakref__, __copy__


Atom
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: Atom
    :exclude-members: __weakref__


Bond
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass :: Bond
    :exclude-members: __weakref__

