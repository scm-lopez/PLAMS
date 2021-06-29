Accessing Old Jobs
==================

The following recipes illustrate how to load data from previously executed jobs:

Binding Native PLAMS Jobs
-------------------------

.. warning::
   The jobs should be loaded with a version of PLAMS that is consistent with the version originally used to run the jobs.


From an existing PLAMS working directory with the contents

::

   OLDDIR/
   ├── OLDJOB1/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB1.dill
   |   ├── OLDJOB1.err
   |   ├── OLDJOB1.in
   |   ├── OLDJOB1.out
   |   ├── OLDJOB1.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── input
   └── logfile

we can bind an instance of the AMSJob class by making use of the `.dill` file.
The AMSJob object in turn contains a results object, which gives access to the data previously calculated.
This can be achieved with the following snippet::

   path       = "OLDDIR/OLDJOB1/OLDJOB1.dill"
   single_JOB = load(path)                                       # AMSJob instance
   if single_JOB.check():
      energy     = single_JOB.results.get_energy()               # load the desired properties
      structure  = single_JOB.results.get_main_molecule()
      propertyX  = single_JOB.results.readrkf('AMSResults', 'DipoleMoment', file='engine')

More often than not, the working directory will include multiple individual subdirectories, each containing individual PLAMS job.

::

   OLDDIR/
   ├── OLDJOB1/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB1.dill
   |   ├── OLDJOB1.err
   |   ├── OLDJOB1.in
   |   ├── OLDJOB1.out
   |   ├── OLDJOB1.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── OLDJOB2/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB2.dill
   |   ├── OLDJOB2.err
   |   ├── OLDJOB2.in
   |   ├── OLDJOB2.out
   |   ├── OLDJOB2.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── OLDJOB3/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB3.dill
   |   ├── OLDJOB3.err
   |   ├── OLDJOB3.in
   |   ├── OLDJOB3.out
   |   ├── OLDJOB3.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── input
   └── logfile

These can be loaded using the `load_all` function and by providing only the path to the top-level directory::

   path       = "OLDDIR"
   all_JOBS   = load_all(path)

Note that `load_all` wraps the `load` function used above and therefore requires existing `.dill` files in each of the loaded subdirectories.
The `load_all` function yields a dictionary with the paths of the `.dill` files as keys and the corresponding job object as values::

   print(all_JOBS)

::

   {'/home/user/OLDDIR/OLDJOB1/OLDJOB1.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baad340b8>,
    '/home/user/OLDDIR/OLDJOB2/OLDJOB2.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baacf24a8>,
    '/home/user/OLDDIR/OLDJOB3/OLDJOB3.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baad06cf8>}

We can now access these AMSJob instances::

   for this_JOB in all_JOBS.values():
      if this_JOB.check():
         energy     = this_JOB.results.get_energy()
         structure  = this_JOB.results.get_main_molecule()
         propertyX  = this_JOB.results.readrkf('AMSResults', 'DipoleMoment', file='engine')


Binding old RKF Files
---------------------
In cases where the `.dill` files are not available any more, it is still possible to load the contents of previously generated `.rkf` files into a PLAMS workflow.
This is done by means of the `KFReader` class::

   path       = "OLDDIR/OLDJOB1/ams.rkf"
   rkf_reader = KFReader(path)
   n_steps    = rkf_reader.read("History", "nEntries")
   energy     = rkf_reader.read("History", "Energy({})".format(n_steps))
   structure  = rkf_reader.read("History", "Coords({})".format(n_steps))

Note that the KFReader class lacks most of the shortcut functions of a proper `AMSResults` object so that the access to the data has to be specified manually.

