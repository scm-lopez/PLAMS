from scm.plams import Settings, JobError, ADFJob, CRSJob, Molecule, ADFResults, ig, CRSResults

__all__ = ['run_crs_adf']


def run_crs_adf(settings_adf, settings_crs,
                solvents, solutes=None,
                return_adfresults=False, **kwargs):
    """A workflow for running COSMO-RS calculations with ADF (*i.e.* DFT) COSMO surface charges.

    The workflow consists of three distinct steps:

        1. Perform gas-phase |ADFJob| calculations on the solvents and solutes (see *settings_adf*).
        2. Perform COSMO |ADFJob| calculations using the .t21 file from step 1. as molecular
           fragment (see *settings_adf*).
           This ensures that zero-point is defined by the gas-phase molecule than the gas-phase
           atomic fragments.
        3. Perform a COSMO-RS calculations with the COSMO surface charges produced in step 2
           (see *settings_crs*).
           This calculation is conducted for all possible solvent/solute pairs,
           assuming solutes have been specified by the user.


    The adf solvation block (*adf_settings.input.solvation*) is soft updated with suitable
    settings for constructing COSMO-RS compatible surface charges
    (see :func:`.add_solvation_block`).
    No user-specified values are overwritten during this process.

    .. admonition:: Examples

        An example value for *settings_adf*:

        .. code:: python

            >>> from scm.plams import Settings

            >>> settings_adf = Settings()
            >>> settings_adf.input.basis.type = 'TZ2P'
            >>> settings_adf.input.xc.gga = 'BP86'
            >>> settings_adf.input.scf.converge = '1.0e-06'

        An example value for *settings_crs* (`activity coefficient`_ calculation):

        .. code:: python

            >>> settings_crs = Settings()
            >>> settings_crs.input.temperature = 298.15
            >>> settings_crs.input.property._h = 'activitycoef'

        And finally the actual calculation with methanol, ethanol and propanol as solvents and
        acetic acid as solute:

        .. code:: python

            >>> solvents = [Molecule('methanol.xyz'), Molecule('ethanol.xyz'), Molecule('propanol.xyz')]
            >>> solutes = Molecule('acetic_acid.xyz')

            >>> crs_dict = run_crs_adf(settings_adf, settings_crs, solvents, solutes)
            >>> print(crs_dict)
            {'CRSJob.methanol.acetic_acid': <scm.plams.interfaces.adfsuite.crs.CRSResults object at 0x7f89b0355668>,
             'CRSJob.ethanol.acetic_acid': <scm.plams.interfaces.adfsuite.crs.CRSResults object at 0x7f89b0355f60>,
             'CRSJob.propanol.acetic_acid': <scm.plams.interfaces.adfsuite.crs.CRSResults object at 0x7f89b0355b00>}

    :type settings_adf: :class:`.Settings`
    :parameter settings_adf:
        A Settings instance with settings for :class:`.ADFJob` (see Examples).

    :type settings_crs: :class:`.Settings`
    :parameter settings_crs:
        A Settings instance with settings for :class:`.CRSJob` (see Examples).

    :type solvents: :class:`.Molecule` or :class:`list` [:class:`.Molecule`]
    :parameter solvents: A Molecule or list of one or more Molecules representing solvents.

    :type solutes: :class:`.Molecule` or :class:`list` [:class:`.Molecule`], optional
    :parameter solutes: An optional Molecule or list of one or more Molecules representing solutes.

    :type return_adfresults: :class:`bool`
    :parameter return_adfresults: If ``True``, return both the solvent and solute ADF results in addition to the final COSMO-RS.

    :parameter \**kwargs, optional:
        Optional keyword arguments that will be passed to all calls of :meth:`.Job.run`.
        For example, one could consider passing a custom jobrunner_ or jobmanager_.

    :returns: A dictionary with the resulting COSMO-RS output.
        The `name` of each :class:`.CRSResults` instance is used as key.
        If ``return_adfresults=True``, return both the COSMO-RS and ADF solvent and solute results.
    :rtype: :class:`dict` or :class:`tuple` [:class:`dict`, :class:`list`, :class:`list`]

    .. _`activity coefficient`: https://www.scm.com/doc/COSMO-RS/Properties.html#activity-coefficients-solvent-and-solute
    .. _jobmanager: ../components/jobmanager.html
    .. _jobrunner: ../components/runners.html

    """  # noqa
    solvents = [solvents] if isinstance(solvents, Molecule) else solvents
    solutes = [solutes] if isinstance(solutes, Molecule) else solutes

    # Validate arguments
    _settings_adf = add_solvation_block(settings_adf)
    _validate_settings_adf(_settings_adf)

    # Decapitalize the "solvation" and "allpoints" keys
    if ig('solvation') in _settings_adf.settings.input:
        _settings_adf.input.solvation = _settings_adf.settings.input.pop(ig('solvation'))
    if ig('allpoints') in _settings_adf.settings.input:
        _settings_adf.input.allpoints = _settings_adf.settings.input.pop(ig('allpoints'))

    # Create the COSMO surfaces for the solute
    solvent_list = [run_adfjob(mol, _settings_adf, **kwargs) for mol in solvents]
    solute_list = [run_adfjob(mol, _settings_adf, **kwargs) for mol in solutes]
    if not solutes:
        solute_list = [None]

    # Start the and return the COSMO-RS job
    crs = {}
    for solvent in solvent_list:
        for solute in solute_list:
            results = run_crsjob(solvent, settings_crs, solute=solute, **kwargs)
            crs[results.job.name] = results

    if not return_adfresults:
        return crs
    else:
        return crs, solvent_list, solute_list


def run_adfjob(mol: Molecule, s: Settings, **kwargs) -> ADFResults:
    """Run an :class:`.ADFJob` on *mol* using the settings provided in *s*."""
    name = 'ADFJob.' + mol.properties.name if 'name' in mol.properties else 'ADFJob.mol'

    # Create the gas-phase molecular fragment
    job1 = ADFJob(molecule=mol, settings=s, name=name+'.gas')
    job1.settings.input.allpoints = ''
    del job1.settings.input.solvation
    results1 = job1.run(**kwargs)
    restart1 = results1['$JN.t21']

    # Construct the ADF COSMO surface
    job2 = ADFJob(molecule=mol, settings=s, depend=[job1], name=name)
    job2.settings.input.allpoints = ''
    job2.settings.input.fragments.gas = restart1
    for at in job2.molecule:
        at.properties.adf.fragment = 'gas'
    return job2.run(**kwargs)


def run_crsjob(solvent: ADFResults, s: Settings, solute: ADFResults = None, **kwargs) -> CRSResults:
    """Run an :class:`.CRSJob` on with *solvent* and, optionally, *solute* using the settings provided in *s*."""
    name = 'CRSJob.' + solvent.job.name.split('.')[1]

    if solute is not None:
        name += '.' + solute.job.name.split('.')[1]
        job = CRSJob(settings=s, depend=[solvent.job, solute.job], name=name)
        set_header(job.settings, solvent['$JN.t21'], solute['$JN.t21'])
    else:
        job = CRSJob(settings=s, depend=[solvent.job], name=name)
        set_header(job.settings, solvent['$JN.t21'])

    return job.run(**kwargs)


def set_header(s: Settings, *values: str) -> None:
    """Assign *value* to the ``["_h"]`` key in *s.input.compound*."""
    s.input.compound = []
    for item in values:
        s.input.compound.append(Settings({'_h': item}))
    s.input.compound[0].frac1 = 1.0  # The first item in *values should be the solvent


def add_solvation_block(adf_settings: Settings) -> None:
    """Add the solvation block to *adf_settings*, returning a copy of the new settings.

    The solvation block (*adf_settings.input.solvation*) is soft updated
    with the following Settings:

    .. code::

        c-mat:  Exact
        charged:        method=Conj
        radii:
            Br:       2.16
            C:        2.0
            Cl:       2.05
            F:        1.72
            H:        1.3
            I:        2.32
            N:        1.83
            O:        1.72
            P:        2.13
            S:        2.16
            Si:       2.48
        scf:    Var All
        solv:   name=CRS cav0=0.0 cav1=0.0
        surf:   Delley

    """
    # Find the solvation key
    solvation = adf_settings.input.find_case('solvation')
    solvation_block = adf_settings.input[solvation]

    # Find all keys for within the solvation block
    keys = ('surf', 'solv', 'charged', 'c-mat', 'scf', 'radii')
    surf, solv, charged, cmat, scf, radii = [solvation_block.find_case(item) for item in keys]

    # Construct the default solvation block
    solvation_block_new = {
        surf: 'Delley',
        solv: 'name=CRS cav0=0.0 cav1=0.0',
        charged: 'method=Conj',
        cmat: 'Exact',
        scf: 'Var All',
        radii: {
            'H': 1.30,
            'C': 2.00,
            'N': 1.83,
            'O': 1.72,
            'F': 1.72,
            'Si': 2.48,
            'P': 2.13,
            'S': 2.16,
            'Cl': 2.05,
            'Br': 2.16,
            'I': 2.32
        }
    }

    # Copy adf_settings and perform a soft update
    ret = adf_settings.copy()
    ret.input[solvation].soft_update(solvation_block_new)
    return ret


def _validate_settings_adf(s: Settings) -> None:
    """Validate the *settings_adf* argument in :func:`run_crs_adf`."""
    solvation = s.input.find_case('solvation')
    solv = s.input[solvation].find_case('solv')

    if solvation not in s.input:
        raise JobError("run_crs_adf: The 'solvation' key is absent"
                       " from settings_adf.input'")

    if solv not in s.input[solvation]:
        raise JobError("run_crs_adf: The 'solv' key is absent"
                       " from settings_adf.input.solvation")

    if 'name=crs' not in s.input[solvation][solv].lower():
        raise JobError("run_crs_adf: The 'name=CRS' value is absent"
                       " from settings_adf.input.solvation.solv")
