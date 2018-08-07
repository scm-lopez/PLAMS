PLAMS
=====

Python Library for Automating Molecular Simulation
------------------------------------------------------

PLAMS is a collection of tools that aims to provide powerful, flexible and easily extendable Python interface to molecular modeling programs. It takes care of input preparation, job execution, file management and output processing as well as helps with building more advanced data workflows.

Usually the daily work of a computational chemist consists of running a number of calculations. Those calculations are done using one or more molecular modeling programs like ADF, BAND, Turbomole or Dirac (we will call such programs *external binaries*). Obtaining results with one of such programs requires a series of steps. First, the subject of the problem (description of a molecular system, set of desired simulation parameters) has to be presented in the format understandable by molecular modeling program and written to an *input file* which is usually a text file. Then the program is executed, it runs and produces *output* which is a collection of text or binary files. That output usually contains more information than is required for a particular problem so data of interest has to be extracted and (possibly) postprocessed. That different computational tools use different input and output formats and are executed differently. In most cases many *single calculations* need to be performed to solve the problem of interest. That requires significant effort to be put into data hygiene to avoid confusing or overwriting input and output files from distinct calculations.

Each of the above steps, apart from actual calculation done by a molecular modeling program, needs to be performed by a human. Preparing and editing text files, creating folders in the filesystem, copying files between them and reading data from output are tedious, repetitive and highly error-prone work. Some users deal with it using automation, usually in form of ad hoc shell scripts. A few programs, like ADF Suite, offer graphical user interface to help with this kind of work, but again, input preparation and output examination, even though assisted with convenient tools, have to be done by a human. Quite often it turns out to be a performance bottleneck to create big  automatic computational workflows, where output data of one calculation is used (usually after some processing) as an input to another calculation, sometimes done with different program on a different machine.

PLAMS was created to solve these problems. It takes responsibility of tiresome and monotonous technical details allowing you to focus on real science and your problem. It lets you do all the things mentioned above (and many more) using simple Python scripts. It gives you a helping hand with automating repetitive or complicated tasks while still leaving you with 100% control over what is really happening with your files, disks and CPUs.


What can be done with PLAMS
----------------------------

The key design principle of PLAMS is *flexibility*. If something (by something we mean: adjusting an input parameter, executing binary with particular options, extracting a value from output etc.) can be done by hand, it can be done with PLAMS. The internal structure of the library was designed in highly modular, object-oriented manner. Thanks to that it takes very little effort to adjust its behavior to one's personal needs or to extend its functionality.


Most important features of PLAMS:
*   preparing, running and examining results of molecular modeling job from within a single Python script
*   convenient automatic file and folder management
*   running jobs in parallel without a need to prepare a special script
*   integration with popular job schedulers (OGE, SLURM, TORQUE)
*   reading and writing molecular coordinates using various formats (`xyz`, `mol`, `mol2`, `pdb`)
*   prevention of multiple runs of the same job
*   easy data transfer between separate runs
*   efficient restarting in case of crash
*   full coverage of all input options and output data in ADF, BAND and DFTB
*   support for Dirac, Orca, Gamess and CP2K (more coming soon)
*   easy extendable for other programs, job schedulers, file formats etc.


Simple example
----------------------------

To provide some real life example: here is a simple PLAMS script which calculates a potential energy curve of a diatomic system:

    #type of atoms (atomic number)
    atom1 = 1
    atom2 = 1

    #interatomic distance values
    dmin = 0.3
    dmax = 1.5
    step = 0.05

    #create a list with interatomic distances
    distances = []
    dist = dmin
    while dist < dmax:
        distances.append(dist)
        dist += step

    #set single calculation parameters (single point, TZ2P/PW91)
    sett = Settings()
    sett.input.basis.type = 'TZ2P'
    sett.input.geometry.sp = True
    sett.input.xc.gga = 'PW91'

    #create a list of jobs
    jobs = []
    for d in distances:
        mol = Molecule()
        mol.add_atom(Atom(atnum=atom1, coords=(0,0,0)))
        mol.add_atom(Atom(atnum=atom2, coords=(d,0,0)))
        job = ADFJob(molecule=mol, settings=sett)
        jobs.append(job)

    #run jobs
    results = [j.run() for j in jobs]

    #extract bond energy from each calculation
    energies = [r.readkf('Energy', 'Bond Energy') for r in results]

    #convert to kcal/mol and print
    energies = [Units.convert(e, 'au', 'kcal/mol') for e in energies]
    print('d[A]    E[kcal/mol]')
    for d,e in zip(distances, energies):
        print('%.2f    %.3f' % (d,e))

When executed, the above script creates uniquely named working folder, then runs 24 independent ADF single point calculations, each in a separate subfolder of the working folder. All files created by each run are saved in the corresponding subfolder for future reference. Finally, the following table describing the potential energy curve of a hydrogen molecule is written to the standard output:

    d[A]    E[kcal/mol]
    0.30    143.301
    0.35    36.533
    0.40    -33.410
    0.45    -79.900
    0.50    -110.823
    0.55    -131.120
    0.60    -143.997
    0.65    -151.598
    0.70    -155.418
    0.75    -156.492
    0.80    -155.572
    0.85    -153.205
    0.90    -149.793
    0.95    -145.635
    1.00    -140.959
    1.05    -135.937
    1.10    -130.699
    1.15    -125.344
    1.20    -119.950
    1.25    -114.576
    1.30    -109.267
    1.35    -104.055
    1.40    -98.967
    1.45    -94.021


Further reading
--------------------

You can find full [PLAMS documentation](https://www.scm.com/doc/plams/index.html) hosted on our website, together with some [tutorials](https://www.scm.com/doc/Tutorials/Scripting/Scripting.html).

You can also build your local copy of the documentation by cloning this repository and executing `doc/build_plams_doc` script.
