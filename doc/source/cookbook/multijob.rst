Dictionary MultiJob
-------------------

Although by default |MultiJob| stores its children in a list, it's often very handy to use a dictionary instead.
A simple example on using a MultiJob to easily organize multiple similar calculations::

    # Read all molecules from the given folder
    molecules = read_molecules('molecules')

    # Calculate, for each molecule, the "total" bond order:
    for molecule in molecules.values():
        molecule.guess_bonds()
        molecule.properties.total_bond_order = sum([bond.order for bond in molecule.bonds])

    # Initialize the common settings:
    settings = Settings()
    settings.input.Basis.Core = 'None'
    settings.input.NumericalQuality = 'Good'
    settings.input.Relativistic = 'Scalar ZORA'
    settings.input.XC.GGA = 'PBE'

    basis = ['SZ', 'DZ', 'DZP', 'TZP', 'TZ2P', 'QZ4P']
    reference_basis = 'QZ4P'

    # Children jobs stored in this MultiJob will be indexed by pairs (molecule_name, basis)
    job = MultiJob(children=dict())

    for bas in basis:
        # Add the basis set to the settings:
        settings.input.Basis.Type = bas
        for name, mol in molecules.items():
            job.children[(name,bas)] = ADFJob(name=name+'_'+bas, molecule=mol, settings=settings)

    job.run()

    # Calculate the average absolute error for each basis set:
    for bas in basis:
        if bas != reference_basis:
            errors = []
            for name, mol in molecules.items():
                ref_energy = job.children[(name,reference_basis)].results.get_energy()
                energy = job.children[(name,bas)].results.get_energy()
                errors.append(abs(energy - ref_energy)/mol.properties.total_bond_order)
            avg_error = Units.convert(sum(errors)/len(errors), 'au', 'kcal/mol')
            print('Basis set: {} Average Absolute Error per bond [kcal/mol]: {}'.format(bas, avg_error))
