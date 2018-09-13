PLAMS cookbook
=========================

MultiJob
--------

A simple example using the MultiJob::

   import os
   from scm.plams import *

   init()

   # Get the molecules:

   molecules_dir = '/absolute/path/molecules'
   molecules = [Molecule(os.path.join(molecules_dir,mol)) for mol in os.listdir(molecules_dir)]

   # Calculate, for each molecule, the "total" bond order:

   total_bond_order = []
   for molecule in molecules:
      molecule.guess_bonds()
      total_bond_order.append(sum([bond.order for bond in molecule.bonds]))

   # Initialize the common settings:

   settings = Settings()

   settings.input.Basis.Core = 'None'
   settings.input.NumericalQuality = 'Good'
   settings.input.Relativistic = 'Scalar ZORA'
   settings.input.XC.GGA = 'PBE'

   basis = ['SZ', 'DZ', 'DZP', 'TZP', 'TZ2P', 'QZ4P']
   reference_basis = 'QZ4P'

   # For each basis set set up a MultiJob and get the bond energies:

   # Let's store the bond energies for different basis sets in a dictionary:
   bond_energies = {}

   for bas in basis:
      # The basis set in the setting object:
      settings.input.Basis.Type = bas

      # Create the MultiJob comprising on  ADFJobs
      job = MultiJob(name=bas, children=[ADFJob(molecule=mol, settings=settings) for mol in molecules])

      # Run the job:
      job.run()

      # Get the bond energies:
      bond_energies[bas] = [child.results.readkf('Energy', 'Bond Energy') for child in job.children]


   # Calculate the average absolute error for each basis set:

   for bas in basis:
      if bas == reference_basis:
         continue
      errs = [abs(ref - en)/n for ref, en, n in zip(bond_energies[reference_basis], bond_energies[bas], total_bond_order)]
      avg_error = sum(errs)/len(errs)
      print('Basis set:', bas, "Average Absolute Error per bond [kcal/mol]:") Units.convert(avg_error, 'au', 'kcal/mol')

   finish()
