COSMO-RS
--------

(*contributed by* `Bas van Beek <https://www.researchgate.net/profile/Bas_van_Beek>`_\)

.. currentmodule:: scm.plams.interfaces.adfsuite.crs

COSMO-RS can be run from PLAMS using the |CRSJOb| class and the corresponding |CRSResults|,
both respectivelly being subclasses of |SCMJob| and |SCMResults|.


Settings
~~~~~~~~

For example, considering the following input block:

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

This input file, for a :math:`sigma`-profile calculation, can be
generated from the following |Settings| instance:

.. code:: python

    >>> print(crs_settings)
    compound:
        _h: /path/to/file.t21
        frac1: 1.0
    property:
        _h: puresigmaprofile
        nprofile: 50
        sigmamax: 0.025
        pure: ''
    temperature: 298.15


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

    >>> compound1, compound2 = Settings(), Settings()
    >>> compound1._h = '/path/to/coumpound1.t21'
    >>> compound1.frac1 = 0.5
    >>> compound2._h = '/path/to/coumpound2.t21'
    >>> compound2.frac1 = 0.5

    >>> s = Settings()
    >>> s.compound = [compound1, compound2]

The resulting code block:

.. code::

    compound /path/to/coumpound1.t21
        frac1 0.5
    end

    compound /path/to/coumpound2.t21
        frac1 0.5
    end

Passing Results object
~~~~~~~~~~~~~~~~~~~~~~

For this reason |Results| instances created by previous jobs can be directly
passed to the ``_h`` key in the ``compound`` block.
The :meth:`CRSJob.run` method will internally attempt to extract any .coskf,
.cos, .crskf or .t21 files from the passed |Results| instance and substitute
itself for the matching path+filename.

.. note::
    :meth:`CRSJob.run` will not check if the extracted file contains any actual
    valid COSMO surface charges.

Example:

.. code:: python

    >>> adf_results = adf_job.run()  # A random ADFJob
    >>> s.compound._h = adf_results
    >>> crs_job = CRSJob(settings=s)

    >>> print(crs_job.settings.compound._h)
    <scm.plams.interfaces.adfsuite.adf.ADFResults object at 0x115444ac8>

    >>> crs_job.run()
    >>> print(crs_job.settings.compound._h)
    /path/to/adf/results.t21

API
~~~

.. autoclass:: CRSJob

.. autoclass:: CRSResults
