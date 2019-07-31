COSMO-RS
--------

(*contributed by* `Bas van Beek <https://www.researchgate.net/profile/Bas_van_Beek>`_\)

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJOb| class and the corresponding |CRSResults|,
both respectivelly being subclasses of |SCMJob| and |SCMResults|.


Settings
~~~~~~~~

For example, considering the following input file for a COSMO-RS
sigma-profile calculation [`1 <https://www.scm.com/doc/COSMO-RS/Analysis.html#sigma-profile>`_]:

.. code::

    compound /path/to/file.t21
        frac1 1.0
    end

    property puresigmaprofile
        nprofile 50
        pure
        sigmamax 0.025
    end

    temperature 298.15

The input file displayed above can be re

.. code:: python

    >>> from scm.plams import Settings, CRSJob

    >>> s = Settings()

    >>> s.input.compound._h = '/path/to/file.t21'
    >>> s.input.compound.frac1 = 1.0
    >>> s.input.property._h = 'puresigmaprofile'
    >>> s.input.property.nprofile = 50
    >>> s.input.property.sigmamax = 0.25
    >>> s.input.property.pure = ''
    >>> s.input.temperature = 298.15

    >>> my_job = CRSJob(settings=s)


Settings with multiple compound
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

More often than not one is interested in the properties of
multi-component mixtures (*e.g.* a dissolved solute).
In such cases one has to pass multiple ``compound`` blocks to
the input file, which is somewhat problematic as Python dictionaries
(including |Settings|) can only contain a set of unique keys.

This problem can be resolved by changing the value of ``compound``
from a |Settings| instance into a list of multiple |Settings| instances.
Each item within this list is expanded into its own ``compound`` block
once :meth:`CRSJob.run` creates the actual input file.

Example |Settings| with two compounds:

.. code:: python

    >>> from scm.plams import Settings, CRSJob

    >>> compound1, compound2, compound3 = Settings(), Settings(), Settings()

    >>> compound1._h = '/path/to/coumpound1.t21'
    >>> compound1.frac1 = 0.33
    >>> compound2._h = '/path/to/coumpound2.t21'
    >>> compound2.frac1 = 0.33
    >>> compound3._h = '/path/to/coumpound3.t21'
    >>> compound3.frac1 = 0.33

    >>> s = Settings()
    >>> s.input.compound = [compound1, compound2, compound3]

    >>> my_job = CRSJob(settings=s)

The job specified above corresponds to the following input file:

.. code::

    compound /path/to/coumpound1.t21
        frac1 0.33
    end

    compound /path/to/coumpound2.t21
        frac1 0.33
    end

    compound /path/to/coumpound3.t21
        frac1 0.33
    end


ADFJob and CRSJob
~~~~~~~~~~~~~~~~~

By default all COSMO-RS solvation energies produced with ADF (*i.e.* DFT) surface charges
are with respect to gas-phase atomic fragments, rather than the gas-phase molecule.
To alleviate this issue a new workflow is presented in the PLAMS cookbook_.
In summary: the construction of the DFT COSMO surface is herein preceded by a restart from
a gas-phase single point. Consequently, the energy of the gas-phase molecule is used as zero point,
rather than the usual gas-phase atomic fragments.

.. _cookbook: ../cookbook/adf_crs.html


API
~~~

.. autoclass:: CRSJob

.. autoclass:: CRSResults
