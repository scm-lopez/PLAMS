ADF (pre-2020 version)
----------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.adf

.. warning::

    This page describes the old interface to the standalone ADF binary.
    Starting from AMS2020, ADF is an AMS engine and should be run using |AMSJob|.
    If you are running AMS2019.3 or older version, you should still use |ADFJob|


ADF can be run from PLAMS using the |ADFJob| class and the corresponding |ADFResults|.
The are subclasses of, respectively, |SCMJob| and |SCMResults|, which gather common pre-AMS logic for all members of the former ADFSuite.


Preparing input
~~~~~~~~~~~~~~~

Input files for ADF consist of blocks and subblocks containg keys and values.
That kind of structure can be easily reflected by |Settings| objects since they are built in a similar way.

The input file is generated based on ``input`` branch of job's |Settings|.
All data present there is translated to input contents.
Nested |Settings| instances define blocks and subblocks, as in the example below::

    myjob = ADFJob(molecule=Molecule('water.xyz'))
    myjob.settings.input.basis.type = 'DZP'
    myjob.settings.input.basis.core = 'None'
    myjob.settings.input.basis.createoutput = 'None'
    myjob.settings.input.scf.iterations = 100
    myjob.settings.input.scf.converge = '1.0e-06 1.0e-06'
    myjob.settings.input.save = 'TAPE13'

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

As you can see, entries present in ``myjob.settings.input.`` are listed in the alphabetical order.
If an entry is a regular key-value pair it is printed in one line (like ``save TAPE13`` above).
If an entry is a nested |Settings| instance it is printed as a block and entries in this instance correspond to contents of a the block.
Both keys and values are kept in their original case.
Strings put as values can contain spaces like ``converge`` above -- the whole string is printed after the key.
That allows to handle lines that need to contain more than one key=value pair.
If you need to put a key without any value, ``True`` or empty string can be given as a value::

    >>> myjob.settings.input.geometry.SP = True
    >>> myjob.settings.input.writefock = ''
    # translates to:
    geometry
      SP
    end

    writefock

If a value of a particualr key is ``False``, that key is omitted.
To produce an empty block simply type::

    >>> myjob.settings.input.geometry  # this is equivalent to myjob.settings.input.geometry = Settings()
    #
    geometry
    end

The algorithm translating |Settings| contents into input file does not check the correctness of the data - it simply takes keys and values from |Settings| instance and puts them in the text file.
Due to that you are not going to be warned if you make a typo, use wrong keyword or improper syntax.
Beware of that.

::

    >>> myjob.settings.input.dog.cat.apple = 'pear'
    #
    dog
      cat
        apple pear
      subend
    end

Some blocks require (or allow) some data to be put in the header line, next to the block name.
Special key ``_h`` is helpful in these situations::

    >>> myjob.settings.input.someblock._h = 'header=very important'
    >>> myjob.settings.input.someblock.key1 = 'value1'
    >>> myjob.settings.input.someblock.key2 = 'value2'
    #
    someblock header=very important
      key1 value1
      key2 value2
    end

The order of blocks within input file and subblocks within a parent block follows |Settings| iteration order which is lexicographical (however, |SCMJob| is smart enough to put blocks like DEFINE or UNITS at the top of the input).
In rare cases you would want to override this order, for example when you supply ATOMS block manually, which can be done when automatic molecule handling is disabled (see below).
That behavior can be achieved by another type of special key::

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

Sometimes one needs to put more instances of the same key within one block, like for example in CONSTRAINTS block in ADF.
It can be done by using list of values instead of a single value::

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

Finally, in some rare cases key and value pair in the input needs to be printed in a form ``key=value`` instead of ``key value``.
When value is a string starting with the equal sign, no space is inserted between key and value::

    >>> myjob.settings.input.block.key = '=value'
    #
    block
      key=value
    end

Sometimes a value of a key in the input file needs to be a path to some file, usually KF file with results of some previous calculation.
Of course such a path can be given explicitly ``newjob.restart = '/home/user/science/plams.12345/oldjob/oldjob.t21'``, but for user's convenience instances of |SCMJob| or |SCMResults| (or directly |KFFile|) can be also used.
Algorithm will detect it and use an absolute path to the main KF file instead::

    >>> myjob.settings.input.restart = oldjob
    >>> myjob.settings.input.fragment.frag1 = fragjob
    #
    restart /home/user/science/plams.12345/oldjob/oldjob.t21
    fragment
      frag1 /home/user/science/fragmentresults/somejob/somejob.t21
    end

|Molecule| instance stored in job's ``molecule`` attribute is automatically processed during the input file preparation and printed in the proper format, depending on the program.
It is possible to disable that and give molecular coordinates explicitly as entries in ``myjob.settings.input.``.
Automatic molecule processing can be turned off by ``myjob.settings.ignore_molecule = True``.



Special atoms in ADF
++++++++++++++++++++

In ADF atomic coordinates in ``atoms`` block can be enriched with some additional information like special names of atoms (for example in case of using different isotopes) or block/fragment membership.
Since usually contents of ``atoms`` block are generated automatically based on the |Molecule| associated with a job, this information needs to be supplied inside the given |Molecule| instance.
Details about every atom can be adjusted separately, by modifying attributes of a particular |Atom| instance according to the following convention:

*   Atomic symbol is generated based on atomic number stored in ``atnum`` attribute of a corresponding |Atom|.
    Atomic number 0 corresponds to the "dummy atom" for which the symbol is empty.
*   If ``atom.properties.ghost`` exists and is ``True``, the above atomic symbol is prefixed with ``Gh.``.
*   If ``atom.properties.name`` exists, its contents are added after the symbol.
    Hence setting ``atnum`` to 0 and adjusting ``name`` allows to put an arbitrary string as the atomic symbol.
*   If ``atom.properties.adf.fragment`` exists, its contents are added after atomic coordinates with ``f=`` prefix.
*   If ``atom.properties.adf.block`` exists, its contents are added after atomic coordinates with ``b=`` prefix.

The following example illustrates the usage of this mechanism::

    >>> mol = Molecule('xyz/Ethanol.xyz')
    >>> mol[1].properties.ghost = True
    >>> mol[2].properties.name = 'D'
    >>> mol[3].properties.ghost = True
    >>> mol[3].properties.name = 'T'
    >>> mol[4].properties.atnum = 0
    >>> mol[4].properties.name = 'J.XYZ'
    >>> mol[5].properties.atnum = 0
    >>> mol[5].properties.name = 'J.ASD'
    >>> mol[5].properties.ghost = True
    >>> mol[6].properties.adf.fragment = 'myfragment'
    >>> mol[7].properties.adf.block = 'block1'
    >>> mol[8].properties.adf.fragment = 'frag'
    >>> mol[8].properties.adf.block = 'block2'
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
~~~~~~~~~~~~~~~~~~~

Runscripts for ADF are very simple.
The only adjustable option (apart from usual ``pre``, ``post``, ``shebang`` and ``stdout_redirect`` which are common for all single jobs) is ``myjob.settings.runscript.nproc``, indicating the number of parallel processes to run ADF with (like with ``-n`` flag or ``NSCM`` environmental variable).





API
~~~

.. autoclass:: ADFResults

Parent abstract classes:

.. currentmodule:: scm.plams.interfaces.adfsuite.scmjob

.. autoclass:: SCMJob(molecule=None, name='plamsjob', settings=None, depend=None)
    :exclude-members: _result_type
.. autoclass:: SCMResults


