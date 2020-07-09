ReaxFF (pre-2019 version)
-------------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.reaxff

.. warning::

    This page describes the old interface to the standalone ReaxFF binary.
    As ReaxFF is now an AMS engine, you probably want to run it using |AMSJob|.

In the few simple examples below you can see how various types of ReaxFF calculations can be performed via PLAMS.

:download:`ReaxFF examples <../_static/reaxff_plams.zip>`


.. autoclass:: ReaxFFJob(molecule=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type


.. autofunction:: load_reaxff_control