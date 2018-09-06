MOPAC
-------------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.mopac


MOPAC (Molecular Orbital PACkage) is a semiempirical quantum chemistry program based on NDDO approximation. More information about MOPAC can be found on its `official website <http://openmopac.net/>`_. MOPAC is also distributed as a part of ADFSuite, where it is integrated with suite's graphical interface.

PLAMS features a basic MOPAC interface defined with classes |MOPACJob| and |MOPACResults|.


Preparing input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Preparing an instance of |MOPACJob| follows general principles for |SingleJob|. Information adjusting input file is stored in ``myjob.settings.input`` branch, whereas a runscript is created based on contents of ``myjob.settings.runscript``. The geometry of your system is supplied in the standard way: with the ``molecule`` attribute.

The input format of MOPAC is simple and straightforward: all the keywords adjusting parameters of your calculation are placed in the first line of the input file. Next two lines are left for user's comments and then geometry of the system follows.

Since blocks and subblocks are not present in MOPAC's input, the ```myjob.settings.input`` branch needs to have a flat structure, just like a regular dictionary, without any nested |Settings| instances. The value of a particular key adjusts the way in which keywords are printed in the first line of the input file:

*   ``myjob.settings.input.keyword = True`` will print ``keyword``
*   ``myjob.settings.input.keyword = value`` will print ``keyword=value`` (with ``value`` being casted to ``str`` if needed)
*   ``myjob.settings.input.keyword = (val1, val2, ...)`` will print ``keyword=(val1,val2,...)`` (when value is a **tuple**)
*   ``myjob.settings.input.keyword = [val1, val2, ...]`` will print ``keyword(val1,val2,...)`` (when value is a **list**)

Moreover, if the keyword ``AUX`` is not supplied by the user, it is automatically inserted in the form ``AUX(0,PRECISSION=9)`` (for compatibility with ADFSuite GUI).

MOPAC allows to freeze each coordinate of each atom separately. This information is extracted from ``mopac_freeze`` key in each atom's properties. If present, this key should contain a string with all axes that you wish to freeze for a particular atom::

    >>> mol = Molecule('system.xyz')
    >>> mol[1].properties.mopac_freeze = 'x'    #freeze x coordinate of atom 1
    >>> mol[2].properties.mopac_freeze = 'yz'   #freeze y and z coordinates of atom 2
    >>> mol[3].properties.mopac_freeze = 'xyz'  #freeze atom 3



API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: MOPACJob(molecule=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type
.. autoclass:: MOPACResults
    :exclude-members: _int2inp
