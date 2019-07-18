COSMO-RS
--------

(*contributed by* `Bas van Beek <https://www.researchgate.net/profile/Bas_van_Beek>`_\)

.. currentmodule:: scm.plams.interfaces.adf.crs

COSMO-RS can be run from PLAMS using the |CRSJOb| class and the corresponding |CRSResults|,
both respectivelly being subclasses of |SCMJob| and |SCMResults|.

Settings
~~~~~~~~

For example, considering t:

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

More often than not one is interested in the properties of
multi-component mixtures (*e.g.* a dissolved solute).
In such cases one has to pass multiple ``compound`` blocks to
the input file, which is somewhat problematic as Python dictionaries
(including |Settings|) can only contain a set of unique keys.

This problem can be resolved by changing the value of ``compound``
from a |Settings| instance into a list of multiple |Settings| instances.
Each item within this list is expanded into its own ``compound`` block
once :meth:`CRSJob.run` creates the actual input file.

For example:

.. code:: python

    >>> compound1 = Settings({'_h': '/path/to/coumpound1.t21', 'frac': 0.25})
    >>> compound2 = Settings({'_h': '/path/to/coumpound2.t21', 'frac': 0.25})
    >>> compound3 = Settings({'_h': '/path/to/coumpound3.t21', 'frac': 0.25})
    >>> compound4 = Settings({'_h': '/path/to/coumpound4.t21', 'frac': 0.25})


API
~~~

.. autoclass:: CRSJob

.. autoclass:: CRSResults
