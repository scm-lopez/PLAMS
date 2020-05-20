from collections import OrderedDict
from ..core.functions import add_to_instance
from ..core.basejob import MultiJob
from ..core.results import Results
from ..core.settings import Settings
from ..mol.molecule import Molecule
from ..mol.atom import Atom
from ..interfaces.adfsuite.ams import AMSJob

__all__ = ['CounterpoiseEnergyJob', 'CounterpoiseEnergyResults']

class CounterpoiseEnergyResults(Results):
    """Results class for counterpoise corrected interaction and binding energies.
    """
    def get_all_energies(self, unit='au'):
        """'A_geo_AB_basis_AB' should be read as "atoms A at the AB geometry with the AB basis set" 
            BSSE is the basis set superposition error
            Eint is the interaction energy
            Ebind is the binding energy (the complex relative to relaxed and isolated A and B)
            _raw means no counterpoise correction
            _cp means counterpoise-corrected
            Edef is the deformation energy
        """
        E = {}
        
        E['AB_geo_AB_basis_AB'] = self.job.children['AB_geo_AB_basis_AB'].results.get_energy(unit=unit) if 'AB_geo_AB_basis_AB' in self.job.children else self.optimized_AB.results.get_energy(unit=unit)
        E['A_geo_A_basis_A'] = self.job.children['A_geo_A_basis_A'].results.get_energy(unit=unit) if 'A_geo_A_basis_A' in self.job.children else self.optimized_A.results.get_energy(unit=unit)
        E['B_geo_B_basis_B'] = self.job.children['B_geo_B_basis_B'].results.get_energy(unit=unit) if 'B_geo_B_basis_B' in self.job.children else self.optimized_B.results.get_energy(unit=unit)
        
            
        E['A_geo_AB_basis_AB'] = self.job.children['A_geo_AB_basis_AB'].results.get_energy(unit=unit)
        E['A_geo_AB_basis_A'] = self.job.children['A_geo_AB_basis_A'].results.get_energy(unit=unit)
        E['B_geo_AB_basis_AB'] = self.job.children['B_geo_AB_basis_AB'].results.get_energy(unit=unit)
        E['B_geo_AB_basis_B'] = self.job.children['B_geo_AB_basis_B'].results.get_energy(unit=unit)


        E['BSSE_A'] = E['A_geo_AB_basis_AB'] - E['A_geo_AB_basis_A']
        E['BSSE_B'] = E['B_geo_AB_basis_AB'] - E['B_geo_AB_basis_B']

        E['BSSE_tot'] = E['BSSE_A'] + E['BSSE_B']

        E['Eint_raw'] = E['AB_geo_AB_basis_AB'] - E['A_geo_AB_basis_A'] - E['B_geo_AB_basis_B']
        E['Eint_cp'] = E['AB_geo_AB_basis_AB'] - E['A_geo_AB_basis_AB'] - E['B_geo_AB_basis_AB']

        E['Edef_A'] = E['A_geo_AB_basis_A'] - E['A_geo_A_basis_A']
        E['Edef_B'] = E['B_geo_AB_basis_B'] - E['B_geo_B_basis_B']

        E['Ebind_raw'] = E['AB_geo_AB_basis_AB'] - E['A_geo_A_basis_A'] - E['B_geo_B_basis_B']
        E['Ebind_cp'] = E['Eint_cp'] + E['Edef_A'] + E['Edef_B']

        return E


class CounterpoiseEnergyJob(MultiJob):
    """A class for calculating the counterpoise-corrected interaction and binding energies.
    """

    _result_type = CounterpoiseEnergyResults

    def __init__(self, ids_state_A, molecule=None, settings_state_AB=None, optimized_AB=None, optimized_A=None, optimized_B=None, settings_state_A=None, settings_state_B=None, **kwargs):
        """
        Specify the complex with either
            - a  molecule and settings_state_AB, or
            - a previously run AMSJob optimized_AB
        You must specify  ids_state_A as a list of one-based atom indices making up the first molecule A. All other atoms are assigned to belong to the second molecule B.
        Optionally set optimized_A and optimized_B to AMSJobs instances of previously run geometry optimizations for the isolated molcules. If they are None, they will be automatically calculated.

        settings_state_A: Setting object for molecule A. Defaults to the same as settings_state_AB. The Task is ignored.
        settings_state_B: Setting object for molecule B. Defaults to the same as settings_state_AB. The Task is ignored.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.optimized_A = optimized_A
        self.optimized_B = optimized_B
        self.optimized_AB = optimized_AB

        # copy the settings so that we wont modify the ones provided as input by the user
        settings_state_AB = settings_state_AB.copy() if settings_state_AB else Settings()
        settings_state_A = settings_state_A.copy() if settings_state_A else settings_state_AB.copy()
        settings_state_B = settings_state_B.copy() if settings_state_B else settings_state_AB.copy()

        # In case the charge key is not specified, explicitly set the value to 0.
        # This is to prevent the charge in molecule.properties.charge (set by get_main_molecule())
        # to be used in case of neutral systems
        for s in [settings_state_AB, settings_state_A, settings_state_B]:
            if not 'charge' in s.input.ams.system:
                s.input.ams.system.charge = 0

        if self.optimized_AB is None:
            self.children['AB_geo_AB_basis_AB'] = AMSJob(molecule, settings=settings_state_AB, name='AB')
            main_child = self.children['AB_geo_AB_basis_AB']
        else:
            main_child = self.optimized_AB
        
        name = 'A_geo_AB_basis_AB'
        s = settings_state_A.copy()
        s.input.ams.Task = 'SinglePoint'
        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            mol = main_child.results.get_main_molecule()
            for i in range(1, len(mol)+1):
                mol[i].properties.ghost = not i in ids_state_A
            self.molecule = mol
        self.children[name] = job

        name = 'A_geo_AB_basis_A'
        s = settings_state_A.copy()
        s.input.ams.Task = 'SinglePoint'
        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            mainmol = main_child.results.get_main_molecule()
            mol = Molecule()
            for i in range(1, len(mainmol)+1):
                if i in ids_state_A:
                    mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
            self.molecule = mol
        self.children[name] = job

        if optimized_A is None:
            name = 'A_geo_A_basis_A'
            s = settings_state_A.copy()
            s.input.ams.Task = 'GeometryOptimization'
            job = AMSJob(settings=s, name=name)
            @add_to_instance(job)
            def prerun(self):
                mainmol = main_child.results.get_main_molecule()
                mol = Molecule()
                for i in range(1, len(mainmol)+1):
                    if i in ids_state_A:
                        mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
                self.molecule = mol
            self.children[name] = job

        name = 'B_geo_AB_basis_AB'
        s = settings_state_B.copy()
        s.input.ams.Task = 'SinglePoint'
        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            mol = main_child.results.get_main_molecule()
            for i in range(1, len(mol)+1):
                mol[i].properties.ghost = i in ids_state_A
            self.molecule = mol
        self.children[name] = job

        name = 'B_geo_AB_basis_B'
        s = settings_state_B.copy()
        s.input.ams.Task = 'SinglePoint'
        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            mainmol = main_child.results.get_main_molecule()
            mol = Molecule()
            for i in range(1, len(mainmol)+1):
                if i not in ids_state_A:
                    mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
            self.molecule = mol
        self.children[name] = job

        if optimized_B is None:
            name = 'B_geo_B_basis_B'
            s = settings_state_B.copy()
            s.input.ams.Task = 'GeometryOptimization'
            job = AMSJob(settings=s, name=name)
            @add_to_instance(job)
            def prerun(self):
                mainmol = main_child.results.get_main_molecule()
                mol = Molecule()
                for i in range(1, len(mainmol)+1):
                    if i not in ids_state_A:
                        mol.add_atom(Atom(atnum=mainmol[i].atnum, coords=mainmol[i].coords))
                self.molecule = mol
            self.children[name] = job


