Disclaimer about AMS2018
~~~~~~~~~~~~~~~~~~~~~~~~~

The 2018 release of Amsterdam Modeling Suite is a very special one. It's the first time the `AMS driver  <../AMS/General.html>`__ unifying the way our software is run is introduced.
However, not all the programs have been yet transformed into AMS engines, which results in a somewhat "hybrid" setup.

Before the 2018 version each of our programs was a separate binary: ADF, BAND, DFTB, ReaxFF, MOPAC, UFF.
The structure of the input file and the form of the results produced were very similar for all of them.
From PLAMS perspective interfaces to all those programs were implemented by introducing a new |SingleJob| subclass called |SCMJob| (and corresponding |Results| subclass called |SCMResults|).
|SCMJob| and |SCMResults| gathered the logic common for all the programs, and were further subclassed into program-specific classes meant to be used in users' scripts: |ADFJob| and |ADFResults|, |BANDJob| and |BANDResults| and so on.

In the 2018 release BAND, DFTB and UFF are no longer separate programs.
They have been turned into engines managed by the common driver called AMS.
As a result of that |BANDJob|, |DFTBJob| and |UFFJob| classes are no longer usable, since they try to run binaries that cannot be found in ``$ADFBIN``.
These classes have been flagged as *legacy interfaces* and will remain in PLAMS for some time for users who still use ADF2017 release (or older).

On the other hand, ADF, ReaxFF, MOPAC and small utilities remain as separate programs (for now).
Some of them (ReaxFF and MOPAC) has been partially integrated with AMS and they can be run both ways: either as a separate program or via the AMS driver.
That means |ADFJob|, |ReaxFFJob|, |MOPACJob| etc. together with the whole |SCMJob| common interface are still default ways of using PLAMS with those programs.

The AMS driver is interfaced from PLAMS with the |AMSJob| class and corresponding |AMSResults|.
These classes are independent from the "old" |SCMJob| class and inherit directly from |SingleJob| and |Results| (although the behavior is often similar to the |SCMJob| family).

.. note::

    How do I run X using PLAMS?

    *   ADF   : |ADFJob|
    *   BAND  : |AMSJob| with BAND engine
    *   DFTB  : |AMSJob| with DFTB engine
    *   ReaxFF: |ReaxFFJob| **or** |AMSJob| with ReaxFF engine (incomplete functionality)
    *   MOPAC : |MOPACJob| **or** |AMSJob| with MOPAC engine
    *   UFF   : |AMSJob| with UFF engine
    *   Densf : |DensfJob|
    *   FCF   : |FCFJob|
