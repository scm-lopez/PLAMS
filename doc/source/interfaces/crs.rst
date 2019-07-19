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


Passing Results object
~~~~~~~~~~~~~~~~~~~~~~

COSMO-RS is somewhat unique among the
For this reason |Results| instances created by previous jobs can be directly
passed to the ``_h`` key in the ``compound`` block.
The :meth:`CRSJob.run` method will internally attempt to extract any .coskf,
.cos, .crskf or .t21 files from the passed |Results| instance and substitute
itself for the matching path+filename.

.. note::
    :meth:`CRSJob.run` will not check if the extracted file contains
    valid COSMO surface charges.

An example is provided below.
In this example the ``adf_results`` represents a |Results| instance
from a previous calculation, like those from a COSMO single point
calculation with |ADFJob|.

.. code:: python

    >>> from scm.plams import Settings, CRSJob

    >>> print(adf_results)  # Results from a random previous job
    <scm.plams.interfaces.adfsuite.adf.ADFResults object at 0x115444ac8>

    >>> s = Settings()
    >>> s.input.compound._h = adf_results
    >>> my_job = CRSJob(settings=s)

    >>> print(my_job.settings.compound._h)
    <scm.plams.interfaces.adfsuite.adf.ADFResults object at 0x115444ac8>

    >>> my_job.run()
    >>> print(my_job.settings.compound._h)
    /path/to/adf/results.t21


ADFJob and CRSJob
~~~~~~~~~~~~~~~~~

.. code:: python

    >>> from scm.plams import ADFJob, CRSJob, Settings

    >>> s_adf = Settings()
    ...

    >>> adf_job1 = ADFJob(settings=s_adf)
    >>> results1 = adf_job1.run()

API
~~~

.. autoclass:: CRSJob

.. autoclass:: CRSResults
