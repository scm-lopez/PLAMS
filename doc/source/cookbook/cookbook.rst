Examples
========

In this chapter we present example PLAMS scripts covering various applications, from very simple tasks like running the same calculation for multiple molecules, to rather advanced dynamic workflows with multiple internal dependencies.

The examples presented in the :ref:`recipes` section are simple job types built using basic PLAMS elements.
They are shipped with PLAMS in the ``recipes`` subpackage and can be directly used in your scripts.
In other words, the code presented there is already included in PLAMS and (unlike examples from two other sections) does not need to be copied to your script.
The source code of ``recipes`` modules is presented here to demonstrate how easy it is to build on top of existing PLAMS elements and create your own fully customized job types.

.. note::

    * PLAMS is an open source library, but in order to run these example scripts you will need a license for the AMS driver and for computational engines used in the examples. Contact license@scm.com for further questions.

    * In order to run these examples the ``ADFBIN`` environment variable should be properly set. You can test this by typing ``$ADFBIN/plams -h`` in a terminal: this should print PLAMS' help message. If this is not the case (e.g. you get 'No such file or directory'), you need to set up the environmental variable ``$ADFBIN`` (see the `Linux Quickstart guide <../../Installation/Linux_Quickstart_Guide.html>`__ for details).

Simple examples
---------------

.. toctree::
 
   WaterOptimization
   ManyJobsInParallel
   ChargeTransferIntegralsADF
   BasisSetBenchmark

Advanced examples
-----------------

.. toctree::

    gammascan

.. _recipes:

Recipes
-------

.. toctree::

    adf_crs
    adffragment
    adfnbo
    numgrad
    numhess
    molecule_gun
    global_minimum
    vibrationASE
