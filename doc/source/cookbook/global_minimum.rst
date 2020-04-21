Global Minimum Search
----------------------

(*contributed by* `Bas van Beek <https://www.researchgate.net/profile/Bas_Beek>`_\)


This module implements a scheme for finding/approaching the conformational global minimum of a |Molecule|. The script accomplishes this by systematically varying all dihedral angles, going through the following steps in the process:

1.    A list of bonds is created, filtering out all bonds are either terminal, part of a ring or have a bond order larger than 1.

2.    The geometry of the molecule is optimized with the first dihedral angle set to 120, 0 and -120 degree; the lowest energy conformer is returned.

3.    This process is repeated in an incremental manner for all valid dihedral angles found in step 1., each step starting from the lowest energy conformer of the previous step.

4.    After all dihedral angles have been exhausted, the final geometry is returned.


Optimizations are possible at various levels of theory, RDKit UFF being a cheap default. Alternatively, the geometry can be optimized at an arbitrary level of theory with the help of the PLAMS |Jobrunner|. Besides the input molecule an additional two arguments, at minimum, are required: A type object of a class derived from |Job| and a dictionary of all keyword arguments that should be passed to aforementioned job (e.g. the job |Settings|).
See below for an exampling using |AMSJob| (DFT level)::

    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.adf.basis.type = 'DZP'
    s.input.adf.XC.GGA = 'PBE'
    s.input.adf.NumericalQuality = 'Good'

    mol_in = Molecule('my_mol.xyz')
    job_kwarg = {'settings': s, 'name': 'my_first_DFTJob'}
    mol_out = global_minimum(mol_in, job_type=AMSJob, **job_kwarg)


Depending on the |Job| class, it may be necessary to manually bind the :meth:`~scm.plams.interfaces.adfsuide.adf.AMSResults.get_energy` and :meth:`~scm.plams.interfaces.adfsuide.adf.AMSResults.get_main_molecule` functions to the jobs matching |Results| class, these two functions being used for reading energies and geometries, respectively.
See below for an exampling using |AMSJob| (DFTB level)::

    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.DFTB.Model = 'DFTB3'
    s.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'

    mol_in = Molecule('my_mol.xyz')
    job_kwarg = {'settings': s, 'name': 'my_first_AMSJob'}
    mol_out = global_minimum(mol_in, job_type=AMSJob, **job_kwarg)


Lastly, by tweaking the job settings including but not excluded to: single points, constrained geometry optimizations, linear transits or a transition state search. The only requirement is that the job yields both an energy and a geometry which can be read with the :meth:`~scm.plams.interfaces.adfsuide.adf.AMSResults.get_energy` and :meth:`~scm.plams.interfaces.adfsuide.adf.AMSResults.get_main_molecule` functions, respectively.
See below for an exampling using |AMSJob| (PBE level). In this example the optimizer searches for a transition state along a reaction coordinate defined by atoms 1 & 2, using the TSRC method, while simultaneously varying all dihedral angles::

    s = Settings()
    s.input.ams.Task = 'TransitionStateSearch'
    s.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = '1 2 1.0'
    s.input.adf.basis.type = 'DZP'
    s.input.adf.XC.GGA = 'PBE'
    s.input.adf.NumericalQuality = 'Good'

    mol_in = Molecule('my_mol.xyz')
    job_kwarg = {'settings': s, 'name': 'my_second_ADF_job'}
    mol_out = global_minimum(mol_in, job_type=AMSJob, **job_kwarg)


The source code:

.. literalinclude:: ../../../recipes/global_minimum.py
