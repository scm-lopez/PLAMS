Charge transfer integrals with ADF
----------------------------------

Calculating charge transfer integrals with ADF requires running a standard two-fragment analysis: two calculations for fragments and one for the full system.
The additional ``TRANSFERINTEGRALS`` key should be used in the full system calculation.
Please refer to the `ADF manual  <../../ADF/Input/Charge_transfer_integrals.html>`_ for more information about charge transfer integrals.

Performing such a 3 step calculation with PLAMS is a straightforward application of :ref:`adf-fragment-recipe` recipe::

    # Add new results extraction method
    @add_to_class(ADFFragmentResults)
    def get_transfer_integrals(self):
        return self.job.full.results._kf.read_section('TransferIntegrals')

    # Common settings for all 3 jobs
    common = Settings()
    common.input.basis.type = 'TZP'
    common.input.basis.core = 'None'
    common.input.basis.createoutput = 'None'
    common.input.xc.gga = 'PW91'
    common.input.symmetry = 'Nosym'
    common.input.noprint = 'Logfile'

    # Specific settings for full system job
    full = Settings()
    full.input.transferintegrals = True

    #Load coordinates from XYZ file and separate into 2 fragments
    mol = Molecule('full_system.xyz')
    mol.guess_bonds()
    fragments = mol.separate()
    if len(fragments) != 2:
        log('ERROR: Molecule {} was split into {} fragments'.format(mol.properties.name, len(fragments)))
        import sys; sys.exit(1)
    else:
        mol1, mol2 = fragments

    # Alternatively one could simply load fragments from separate xyz files:
    # mol1 = Molecule('fragment1.xyz')
    # mol2 = Molecule('fragment2.xyz')

    job = ADFFragmentJob(name='ADFTI',fragment1=mol1,fragment2=mol2,settings=common,full_settings=full)
    results = job.run()

    #TI is a dictionary with the whole TransferIntegrals section from TAPE21
    TI = results.get_transfer_integrals()
    for key in TI:
        print(key, TI[key])

In the code above you can see the usage of |binding_decorators| to enrich ``ADFFragmentResults`` with a method extracting transfer integrals.
Another feature worth noting is automatic separation of a molecular system into two fragments.
