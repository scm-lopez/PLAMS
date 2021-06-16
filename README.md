PLAMS
=====

Python Library for Automating Molecular Simulation
------------------------------------------------------

PLAMS is a collection of tools that aims to provide powerful, flexible and easily extendable Python interface to molecular modeling programs. It takes care of input preparation, job execution, file management and output processing as well as helps with building more advanced data workflows.

Usually the daily work of a computational chemist consists of running a number of calculations. Those calculations are done using one or more molecular modeling programs like ADF, BAND, Turbomole or Dirac (we will call such programs *external binaries*). Obtaining results with one of such programs requires a series of steps. First, the subject of the problem (description of a molecular system, set of desired simulation parameters) has to be presented in the format understandable by molecular modeling program and written to an *input file* which is usually a text file. Then the program is executed, it runs and produces *output* which is a collection of text or binary files. That output usually contains more information than is required for a particular problem so data of interest has to be extracted and (possibly) postprocessed. That different computational tools use different input and output formats and are executed differently. In most cases many *single calculations* need to be performed to solve the problem of interest. That requires significant effort to be put into data hygiene to avoid confusing or overwriting input and output files from distinct calculations.

Each of the above steps, apart from actual calculation done by a molecular modeling program, needs to be performed by a human. Preparing and editing text files, creating folders in the filesystem, copying files between them and reading data from output are tedious, repetitive and highly error-prone work. Some users deal with it using automation, usually in form of ad hoc shell scripts. A few programs, like the Amsterdam Modeling Suite, offer graphical user interface to help with this kind of work, but again, input preparation and output examination, even though assisted with convenient tools, have to be done by a human. Quite often it turns out to be a performance bottleneck to create big  automatic computational workflows, where output data of one calculation is used (usually after some processing) as an input to another calculation, sometimes done with different program on a different machine.

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
*   full coverage of all input options and output data for programs of the Amsterdam Modeling Suite (ADF, BAND, DFTB, MOPAC, ReaxFF, ...)
*   support for Dirac, Orca, Gamess and CP2K (more coming soon)
*   easy extendable for other programs, job schedulers, file formats etc.


Simple example
----------------------------

To provide some real life example, here is a simple PLAMS script performing the following work-flow:

- generate xyz geometries from SMILES strings
- optimize the structures using the fast semi-empirical method GFN1-xTB
- compute the bonding energy for the optimized structures using the DFT engine ADF

```python
# A dictionary with molecule names and smiles strings:
mol_smiles = {'Methane'  : 'C',
              'Ethane'   : 'C-C',
              'Ethylene' : 'C=C',
              'Acetylene': 'C#C'}

# Settings for the semi-empirical GFN1-xTB geometry optimization:
go_xtb_settings = Settings()
go_xtb_settings.input.ams.Task = 'GeometryOptimization'
go_xtb_settings.input.dftb.Model = 'GFN1-xTB'

# Settings for the single point DFT calculation with ADF:
sp_adf_settings = Settings()
sp_adf_settings.input.ams.Task = 'SinglePoint'
sp_adf_settings.input.adf.basis.type = 'DZP'
sp_adf_settings.input.adf.XC.GGA = 'PBE'

for name, smiles in mol_smiles.items():
    # Generate an xyz structure from the smiles string:
    mol = from_smiles(smiles)

    # Run the geometry optimization with GFN1-xTB:
    go_job = AMSJob(molecule=mol, settings=go_xtb_settings, name=f'go_{name}')
    go_job.run()
    optimized_mol = go_job.results.get_main_molecule()

    # Run the single point energy calcualtion with DFT (ADF:
    sp_job = AMSJob(molecule=optimized_mol, settings=sp_adf_settings, name=f'sp_{name}')
    sp_job.run()
    bonding_energy = sp_job.results.get_energy(unit='kcal/mol')

    # Print the results:
    print(f"Energy for {name}: {bonding_energy} [kcal/mol]")
```

When executed, the above script creates uniquely named working folder, then runs a series of calculations, each in a separate subfolder of the working folder. All files created by each run are saved in the corresponding subfolder for future reference. 


Further reading
--------------------

You can find full [PLAMS documentation](https://www.scm.com/doc/plams/index.html) hosted on our website, together with some [tutorials](https://www.scm.com/doc/Tutorials/Scripting/Scripting.html).

You can also build your local copy of the documentation by cloning this repository and executing `doc/build_plams_doc` script.
