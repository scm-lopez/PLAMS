ADF Suite
-------------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.scmjob

PLAMS offers interfaces to three main binaries of the ADF Suite: ADF, BAND and DFTB as well as some other small utility binaries like DENSF of FCF. All possible input keywords and options are covered, as well as extraction of arbitrary data from binary files (called KF files) produced by these programs.

.. _adf-band-dftb:

ADF, BAND and DFTB
~~~~~~~~~~~~~~~~~~~~~~~~~

ADF, BAND and DFTB are of course very different programs, but from our perspective they are rather similar. Their input files follow a common structure of blocks and subblocks. They store results as binary files in KF format and and they print human-readable summary of calculation to the standard output. They also share command line arguments, error messages etc. Thanks to that Python code responsible for creating, running and examining jobs with ADF, BAND and DFTB jobs overlaps a lot and can be grouped together in abstract classes. |SCMJob| and |SCMResults| are subclasses of, respectively, |SingleJob| and |Results| and serve as bases for concrete classes: |ADFJob|, |BANDJob|, |DFTBJob|, |ADFResults|, |BANDResults| and |DFTBResults|. Code contained in these concrete classes describes small technical differences and is used only internally, so they are omitted in the API specification below. From user perspective they all follow the common interface defined by |SCMJob| and |SCMResults|. That means in your scripts you would create instances of |ADFJob|, |BANDJob| or |DFTBJob|, but methods that you can use with them (and their corresponding results) can be taken from |SCMJob| and |SCMResults|.



Preparing input
+++++++++++++++

Although input files for ADF, BAND and DFTB use different sets of keywords, they all have the same logical structure -- they consist of blocks and subblocks containg keys and values. That kind of structure can be easily reflected by |Settings| objects since they are built in a similar way.

The input file is generated based on ``input`` branch of job's |Settings|. All data present there is translated to input contents. Nested |Settings| instances define blocks and subblocks, as in the example below::

    >>> myjob = ADFJob(molecule=Molecule('water.xyz'))
    >>> myjob.settings.input.basis.type = 'DZP'
    >>> myjob.settings.input.basis.core = 'None'
    >>> myjob.settings.input.basis.createoutput = 'None'
    >>> myjob.settings.input.scf.iterations = 100
    >>> myjob.settings.input.scf.converge = '1.0e-06 1.0e-06'
    >>> myjob.settings.input.save = 'TAPE13'

Input file created during execution of ``myjob`` looks like::

    atoms
        #coordinates from water.xyz
    end

    basis
      createoutput None
      core None
      type DZP
    end

    save TAPE13

    scf
      converge 1.0e-06 1.0e-06
      iterations 100
    end

As you can see, entries present in ``myjob.settings.input.`` are listed in the alphabetical order. If an entry is a regular key-value pair it is printed in one line (like ``save TAPE13`` above). If an entry is a nested |Settings| instance it is printed as a block and entries in this instance correspond to contents of a the block. Both keys and values are kept in their original case. Strings put as values can contain spaces like ``converge`` above -- the whole string is printed after the key. That allows to handle lines that need to contain more than one key=value pair. If you need to put a key without any value, ``True`` or empty string can be given as a value::

    >>> myjob.settings.input.geometry.SP = True
    >>> myjob.settings.input.writefock = ''
    # translates to:
    geometry
      SP
    end

    writefock

If a value of a particualr key is ``False``, that key is omitted. To produce an empty block simply type::

    >>> myjob.settings.input.geometry  # this is equivalent to myjob.settings.input.geometry = Settings()
    #
    geometry
    end

The algorithm translating |Settings| contents into input file does not check the correctness of the data - it simply takes keys and values from |Settings| instance and puts them in the text file. Due to that you are not going to be warned if you make a typo, use wrong keyword or improper syntax. Beware of that.

::

    >>> myjob.settings.input.dog.cat.apple = 'pear'
    #
    dog
      cat
        apple pear
      subend
    end

Some blocks require (or allow) some data to be put in the header line, next to the block name. Special key ``_h`` is helpful in these situations::

    >>> myjob.settings.input.someblock._h = 'header=very important'
    >>> myjob.settings.input.someblock.key1 = 'value1'
    >>> myjob.settings.input.someblock.key2 = 'value2'
    #
    someblock header=very important
      key1 value1
      key2 value2
    end

The order of blocks within input file and subblocks within a parent block follows |Settings| iteration order which is lexicographical (however, |SCMJob| is smart enough to put blocks like DEFINE or UNITS at the top of the input). In rare cases you would want to override this order, for example when you supply ATOMS block manually, which can be done when automatic molecule handling is disabled (see below). That behavior can be achieved by another type of special key::

    >>> myjob.settings.input.block._1 = 'entire line that has to be the first line of block'
    >>> myjob.settings.input.block._2 = 'second line'
    >>> myjob.settings.input.block._4 = 'I will not be printed'
    >>> myjob.settings.input.block.key1 = 'value1'
    >>> myjob.settings.input.block.key2 = 'value2'
    #
    block
      entire line that has to be the first line of block
      second line
      key1 value1
      key2 value2
    end

Sometimes one needs to put more instances of the same key within one block, like for example in CONSTRAINTS block in ADF. It can be done by using list of values instead of a single value::

    >>> myjob.settings.input.constraints.atom = [1,5,4]
    >>> myjob.settings.input.constraints.block = ['ligand', 'residue']
    #
    constraints
      atom 1
      atom 5
      atom 4
      block ligand
      block residue
    end

Finally, in some rare cases key and value pair in the input needs to be printed in a form ``key=value`` instead of ``key value``. When value is a string starting with the equal sign, no space is inserted between key and value::

    >>> myjob.settings.input.block.key = '=value'
    #
    block
      key=value
    end

Sometimes a value of a key in the input file needs to be a path to some file, usually KF file with results of some previous calculation. Of course such a path can be given explicitly ``newjob.restart = '/home/user/science/plams.12345/oldjob/oldjob.t21'``, but for user's convenience instances of |SCMJob| or |SCMResults| (or directly |KFFile|) can be also used. Algorithm will detect it and use an absolute path to the main KF file instead::

    >>> myjob.settings.input.restart = oldjob
    >>> myjob.settings.input.fragment.frag1 = fragjob
    #
    restart /home/user/science/plams.12345/oldjob/oldjob.t21
    fragment
      frag1 /home/user/science/fragmentresults/somejob/somejob.t21
    end

|Molecule| instance stored in job's ``molecule`` attribute is automatically processed during the input file preparation and printed in the proper format, depending on the program. It is possible to disable that and give molecular coordinates explicitly as entries in ``myjob.settings.input.``. Automatic molecule processing can be turned off by ``myjob.settings.ignore_molecule = True``.



Special atoms in ADF
++++++++++++++++++++

In ADF atomic coordinates in ``atoms`` block can be enriched with some additional information like special names of atoms (for example in case of using different isotopes) or block/fragment membership. Since usually contents of ``atoms`` block are generated automatically based on the |Molecule| associated with a job, this information needs to be supplied inside the given |Molecule| instance. Details about every atom can be adjusted separately, by modifying attributes of a particular |Atom| instance according to the following convention:

*   Atomic symbol is generated based on atomic number stored in ``atnum`` attribute of a corresponding |Atom|. Atomic number 0 corresponds to the "dummy atom" for which the symbol is empty.
*   If an attribute ``ghost`` of an |Atom| is ``True``, the above atomic symbol is prefixed with ``Gh.``.
*   If an |Atom| has an attribute ``name`` its contents are added after the symbol. Hence setting ``atnum`` to 0 and adjusting ``name`` allows to put an arbitrary string as the atomic symbol.
*   If an |Atom| has an attribute ``fragment`` its contents are added after atomic coordinates with ``f=`` prefix.
*   If an |Atom| has an attribute ``block`` its contents are added after atomic coordinates with ``b=`` prefix.

The following example illustrates the usage of this mechanism::

    >>> mol = Molecule('xyz/Ethanol.xyz')
    >>> mol[0].ghost = True
    >>> mol[1].name = 'D'
    >>> mol[2].ghost = True
    >>> mol[2].name = 'T'
    >>> mol[3].atnum = 0
    >>> mol[3].name = 'J.XYZ'
    >>> mol[4].atnum = 0
    >>> mol[4].name = 'J.ASD'
    >>> mol[4].ghost = True
    >>> mol[5].fragment = 'myfragment'
    >>> mol[6].block = 'block1'
    >>> mol[7].fragment = 'frag'
    >>> mol[7].block = 'block2'
    >>> myjob = ADFJob(molecule=mol)
    #
    atoms
          1      Gh.C       0.01247       0.02254       1.08262
          2       C.D      -0.00894      -0.01624      -0.43421
          3    Gh.H.T      -0.49334       0.93505       1.44716
          4     J.XYZ       1.05522       0.04512       1.44808
          5  Gh.J.ASD      -0.64695      -1.12346       2.54219
          6         H       0.50112      -0.91640      -0.80440 f=myfragment
          7         H       0.49999       0.86726      -0.84481 b=block1
          8         H      -1.04310      -0.02739      -0.80544 f=frag b=block2
          9         O      -0.66442      -1.15471       1.56909
    end





Preparing runscript
+++++++++++++++++++

Runscripts for ADF, BAND and DFTB are very simple - they are just single execution of one of the binaries with proper standard input and output handling. The number of parallel processes (``-n`` parameter) can be adjusted with ``myjob.settings.runscript.nproc``.



Results extraction
++++++++++++++++++

All three programs print results to the standard output. The output file can be examined with standard text processing tools (:meth:`~scm.plams.core.results.Results.grep_output` and :meth:`~scm.plams.core.results.Results.awk_output`). Besides that all calculation details are saved in the binary file in KF format. This file is called ``TAPE21`` for ADF, ``RUNKF`` for BAND and ``dftb.rkf`` for DFTB. PLAMS renames those files to, respectively ``[jobname].t21``, ``[jobname].runkf`` and ``[jobname].rkf``. Data stored in those files can be accessed using additional methods defined in |SCMResults| class.


API
+++

.. autoclass:: SCMJob(molecule=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type
.. autoclass:: SCMResults



Other tools: Densf, FCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from main computational programs mentioned above, ADFSuite offers a range of small utility tools that can be used to obtain more specific results. Those tools usually base on the prior run of one of the main programs and need the KF file produced by them as a part of the input.

From the functional point of view these tools are very similar to ADF, BAND and DFTB. Their results are stored in KF files and their input files follow the same structure of blocks, keys and values. Because of that the same classes (|SCMJob| and |SCMResults|) are used as bases and hence preparation, running and results extraction for utility tools follow the rules described above, in :ref:`adf-band-dftb`

The main difference is that usually utility jobs don't need molecular coordinates as part of the input (they extract this information from previous calculation's KF file). So no |Molecule| instance is needed and the ``molecule`` attribute of the job object is simply ignored. Because of that :meth:`~SCMResults.get_molecule` method does not work with :class:`FCFResults`, :class:`DensfResults` etc.

Below you can find the list of dedicated job classes that are currently available. Details about input specification for those jobs can be found in corresponding part of ADF Suite documentation.


.. currentmodule:: scm.plams.interfaces.adfsuite.fcf

.. autoclass:: FCFJob(inputjob1=None, inputjob2=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type


.. currentmodule:: scm.plams.interfaces.adfsuite.densf

.. autoclass:: DensfJob(inputjob=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type

.. _kf_files:

KF files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: scm.plams.tools.kftools

KF is the main format for storing binary data used in all ADFSuite programs. PLAMS offers an easy and efficient way of accessing the data stored in existing KF files, as well as modifying and creating them.

KFFile
++++++

.. autoclass:: KFFile
    :exclude-members: __weakref__

KFReader
++++++++

.. autoclass:: KFReader
    :exclude-members: __weakref__
