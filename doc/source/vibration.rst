Vibration
-------------------------

.. currentmodule:: scm.plams.recipes.vibration

If the ASE module is present, it can be used to calculate numerical frequencies in parallel.
See the `ASE documentation <https://wiki.fysik.dtu.dk/ase/ase/vibrations/vibrations.html>`_ for a detailed overview of its capabilities. The |VibrationsJob| class provides an easy to use interface to the ``ase.vibrations.Vibrations`` class, allowing for parallel calculations of single-point gradients and therefore speeding up any Frequency calculation by a factor equal to the atomic displacements beeing calculated.


Example
~~~~~~~~~~~~~~~~~~~~~~~~~
::

            settings = Settings()
            settings.input.symmetry = 'nosym'
            settings.input.basis.type = 'DZ'
            settings.input.basis.core = 'None'
            settings.input.xc.gga = 'PBE'
            settings.input.basis.createoutput = 'None'
            settings.input.noprint = 'logfile'
            settings.input.gradient = True
            #Use previously converged coefficients to speed up convergence - saves a lot of CPU time!!
            settings.input.restart._h = os.path.join(os.getcwd(),'restart.t21')+' &'
            settings.input.restart.nogeo = True

            mol = Molecule('mol.xyz') # read Molecule from mol.xyz

            vibJob = VibrationsJob(molecule=mol, settings=settings, aseVibOpt={ 'indices': [1], 'nfree': 4 })

            vibJob.run()

            vib = vibJob.results.get_ASEVib()
            vib.summary(method='Frederiksen')



API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: VibrationsJob(molecule, settings, jobType=ADFJob, get_gradients='get_gradients', reorder='inputOrder', aseVibOpt={})
    :exclude-members: _result_type
.. autoclass:: VibrationsResults
