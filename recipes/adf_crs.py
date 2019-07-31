from scm.plams import Settings, JobError, ADFJob, CRSJob

__all__ = ['run_crs_adf']


def run_crs_adf(mol, settings_adf, settings_crs, **kwargs):
    """A workflow for running COSMO-RS calculations with ADF (*i.e.* DFT) COSMO surface charges.

    The workflow consists of three distinct steps:

        1. Perform a gas-phase |ADFJob| calculation (see *settings_adf*).
        2. Perform a COSMO |ADFJob| calculation using the .t21 file from step 1. as molecular
           fragment (see *settings_adf*).
           This ensures that zero-point is defined by the gas-phase molecule than the gas-phase
           atomic fragments.
        3. Perform a COSMO-RS calculation with the COSMO surface charges produced in step 2
           (see *settings_crs*).

    If the calculation involves multiple components (*e.g.* a solvent and solute),
    than the compound header (``"_h"``) of the solute will be recognized by its value: ``None``
    (see the *settings_crs* example below).
    A :exc:`.JobError` is raised if no such value is specified.

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

            >>> from scm.plams import Settings

            >>> solvent, solute = Settings(), Settings()
            >>> solvent._h = '/path/to/solvent.t21'
            >>> solvent.frac1 = 1.0

            # run_crs_adf() will automatically replace None with the to-be created solute .t21 file
            >>> solute._h = None

            >>> settings_crs = Settings()
            >>> settings_crs.input.compound = [solvent, solute]
            >>> settings_crs.input.temperature = 298.15
            >>> settings_crs.input.property._h = 'activitycoef'

    :type mol: :class:`.Molecule`
    :parameter mol: A PLAMS Molecule.

    :type settings_adf: :class:`.Settings`
    :parameter settings_adf:
        A Settings instance with settings for :class:`.ADFJob` (see Examples).

    :type settings_crs: :class:`.Settings`
    :parameter settings_crs:
        A Settings instance with settings for :class:`.CRSJob` (see Examples).

    :parameter \**kwargs:
        Optional keyword arguments that will be passed to all calls of :meth:`.Job.run`.
        For example, one could consider passing a custom jobrunner_ or jobmanager_.

    :returns: The resulting COSMO-RS output.
    :rtype: :class:`.CRSResults`

    .. _`activity coefficient`: https://www.scm.com/doc/COSMO-RS/Properties.html#activity-coefficients-solvent-and-solute
    .. _jobmanager: ../components/jobmanager.html
    .. _jobrunner: ../components/runners.html

    """  # noqa
    # Validate arguments
    _settings_adf = add_solvation_block(settings_adf)
    _validate_settings_adf(_settings_adf)
    _validate_settings_crs(settings_crs)

    # Identify the "solvation" and "allpoints" keys
    solvation = _settings_adf.input.find_case('solvation')
    allpoints = _settings_adf.input.find_case('allpoints')

    # Create the gas-phase molecular fragment
    adf_job1 = ADFJob(molecule=mol, settings=_settings_adf)
    adf_job1.settings.input[allpoints] = ''
    del adf_job1.settings.input[solvation]
    adf_results1 = adf_job1.run(**kwargs)
    adf_restart1 = adf_results1['$JN.t21']

    # Construct the ADF COSMO surface
    adf_job2 = ADFJob(molecule=mol, settings=_settings_adf, depend=[adf_job1])
    adf_job2.settings.input[allpoints] = ''
    adf_job2.settings.input.fragments.gas = adf_restart1
    for at in adf_job2.molecule:
        at.properties.adf.fragment = 'gas'
    adf_results2 = adf_job2.run(**kwargs)
    adf_restart2 = adf_results2['$JN.t21']

    # Start the and return the COSMO-RS job
    crs_job = CRSJob(settings=settings_crs, depend=[adf_job1, adf_job2])
    set_header(crs_job.settings, adf_restart2)
    ret = crs_job.run(**kwargs)
    ret.wait()
    return ret


def set_header(s: Settings, value: str) -> None:
    """Assign *value* to the ``["_h"]`` key in *s.input.compound*."""
    compound = s.input.find_case('compound')
    compound_block = s.input[compound]

    if isinstance(compound_block, Settings):
        compound_block._h = value
    else:
        # Find the first element where the _h key is None
        for item in compound_block:
            if item._h is None:
                item._h = value
                break


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


def _validate_settings_crs(s: Settings) -> None:
    """Validate the *settings_crs* argument in :func:`run_crs_adf`."""
    compound = s.input.find_case('compound')
    value = s.input[compound]

    if compound not in s.input:
        raise JobError("run_crs_adf: The 'compound' key is absent from settings_crs.input")

    if not isinstance(value, (list, Settings)):
        class_name = repr(value.__class__.__name__)
        raise JobError("run_crs_adf: settings_crs.input.compound expects an instance "
                       "of 'list' or 'Settings'; observed type: " + class_name)

    if isinstance(value, list):
        header_list = [item._h for item in value]
        if header_list.count(None) != 1:
            raise JobError("run_crs_adf: settings_crs.input.compound was passed as a 'list'; "
                           "a single 'None' was expected in the '_h' key of one of its values")
