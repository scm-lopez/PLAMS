from scm.plams import Settings, JobError, ADFJob, CRSJob


def run_crs_adf(mol, settings_adf, settings_crs, **kwargs):
    """A workflow for running COSMO-RS calculations using ADF COSMO surface charges.

    The workflow consists of three distinct steps:

        1. Perform a gas-phase |ADFJob| calculation (see *settings_adf*).
        2. Perform a COSMO |ADFJob| calculation using the .t21 file from step 1. as molecular
           fragment (see *settings_adf*).
           This ensures that the total bonding energy produced in step 1. is used as zero-point,
           *i.e.* the energy reported in step 2. is the actual solvation energy.
        3. Perform a COSMO-RS calculation using the surface charges produced in step 2
           (see *settings_crs*).

    :type mol: :class:`.Molecule`
    :parameter mol: A PLAMS Molecule.

    :type settings_adf: :class:`.Settings`
    :parameter settings_adf:
        A Settings instance with settings for :class:`.ADFJob`.

    :type settings_crs: :class:`.Settings`
    :parameter settings_crs:
        A Settings instance with settings for :class:`.CRSJob`.

    :parameter \**kwargs:
        Optional keyword arguments that will be passed to all calls of the :meth:`.Job.run` method.
        For example, one could pass a custom *jobrunner* or *jobmanager*.

    :returns: A :class:`.CRSResults` instance.

    """  # noqa
    _validate_settings_adf(settings_adf)
    _validate_settings_crs(settings_crs)
    solvation = settings_adf.input.find_case('solvation')
    allpoints = settings_adf.input.find_case('allpoints')

    # Create the gas-phase molecular fragment
    adf_job1 = ADFJob(mol, settings=settings_adf)
    adf_job1.settings.input[allpoints] = ''
    del adf_job1.settings.input[solvation]
    adf_results1 = adf_job1.run(**kwargs)
    adf_restart1 = adf_results1['$JN.t21']

    # Construct the ADF COSMO surface
    adf_job2 = ADFJob(mol, settings=settings_adf, depend=adf_job1)
    adf_job2.settings.input[allpoints] = ''
    adf_job2.settings.input.fragments.gas = adf_restart1
    for at in adf_job2.molecule:
        at.properties.adf.fragment = 'gas'
    adf_results2 = adf_job2.run(**kwargs)
    adf_restart2 = adf_results2['$JN.t21']

    # Start the and return the COSMO-RS job
    crs_job = CRSJob(settings=settings_crs, depend=adf_job2)
    set_header(crs_job.settings, adf_restart2)
    ret = crs_job.run(**kwargs)
    ret.wait()
    return ret


def set_header(s, value):
    """Assign *value* to the ``["_h"]`` key in *s.input.compound*."""
    compound = s.input.find_case('compound')

    if isinstance(s, Settings):
        s.input[compound]._h = value
    else:
        # Find the first element where the _h key is None
        for item in s.input[compound]:
            if item._h is None:
                item._h = value
                break


def _validate_settings_adf(s):
    """Validate the *settings_adf* argument in :func:`run_crs_adf`."""
    solvation = s.input.find_case('solvation')
    solv = s.input[solvation].find_case('solv')

    if solvation not in s.input:
        raise JobError("run_crs_adf: The 'solvation' key is absent"
                       " from settings_adf.input'")

    if solv not in s.input[solvation]:
        raise JobError("run_crs_adf: The 'solv' key is absent"
                       " from settings_adf.input.solvation'")

    if 'name=crs' not in s.input[solvation][solv].lower():
        raise JobError("run_crs_adf: The 'name=CRS' value is absent"
                       " from settings_adf.input.solvation.solv")


def _validate_settings_crs(s):
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
