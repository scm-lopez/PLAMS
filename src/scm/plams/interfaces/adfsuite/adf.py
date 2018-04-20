from .scmjob import SCMJob, SCMResults
from ...core.errors import ResultsError
from ...tools.units import Units
from ...tools.periodic_table import PT


__all__ = ['ADFJob', 'ADFResults']


class ADFResults(SCMResults):
    _kfext = '.t21'
    _rename_map = {'TAPE{}'.format(i) : '$JN.t{}'.format(i) for i in range(10,100)}


    def get_properties(self):
        """get_properties()
        Return a dictionary with all the entries from ``Properties`` section in the main KF file (``$JN.t21``).
        """
        ret = {}
        for sec,var in self._kf:
            if sec == 'Properties':
                ret[var] = self.readkf(sec,var)
        return ret

    def get_main_molecule(self):
        """get_main_molecule()
        Return a |Molecule| instance based on the ``Geometry`` section in the main KF file (``$JN.t21``).

        For runs with multiple geometries (geometry optimization, transition state search, intrinsic reaction coordinate) this is the **final** geometry. In such a case, to access the initial (or any intermediate) coordinates please use :meth:`get_input_molecule` or extract coordinates from section ``History``, variables ``xyz 1``, ``xyz 2`` and so on. Mind the fact that all coordinates written by ADF to ``History`` section are in bohr and internal atom order::

            mol = results.get_molecule(section='History', variable='xyz 1', unit='bohr', internal=True)
        """
        return self.get_molecule(section='Geometry', variable='xyz InputOrder', unit='bohr')


    def get_input_molecule(self):
        """get_input_molecule()
        Return a |Molecule| instance with initial coordinates.

        All data used by this method is taken from ``$JN.t21`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        if ('History', 'nr of geometries') in self._kf:
            return self.get_molecule(section='History', variable='xyz 1', unit='bohr', internal=True)
        return self.get_main_molecule()


    def get_energy(self, unit='au'):
        """get_energy(unit='au')
        Return final bond energy, expressed in *unit*.
        """
        return self._get_single_value('Energy', 'Bond Energy', unit)


    def get_dipole_vector(self, unit='au'):
        """get_dipole_vector(unit='au')
        Return the dipole vector, expressed in *unit*.
        """
        prop = self.get_properties()
        if 'Dipole' in prop:
            return Units.convert(prop['Dipole'], 'au', unit)
        raise ResultsError("'Dipole' not present in 'Properties' section of {}".format(self._kfpath()))


    def get_ordered_gradients(self, eUnit='a.u.', lUnit='bohr'):
        """get_ordered_gradients
        Return cartesian gradients in order of the input molecule. Returns a numpy array
        with shape (nAtoms,3) in the units given (standard a.u./bohr).
        """
        from numpy import array as npAr
        gradients = self.get_gradients(eUnit=eUnit, lUnit=lUnit)
        return npAr(self.inputOrder(gradients))

    def get_gradients(self, eUnit='a.u.', lUnit='bohr'):
        """get_gradients
        Returns the cartesian gradients from the 'Gradients_CART' field of the 'GeoOpt' Section in the kf-file
        as a numpy array with shape (nAtoms,3) in the units given (standard a.u./bohr).
        Note: The values from the KF-File are multiplied with -1.0, reason unknown.
        """
        from numpy import array as npAr
        from numpy import reshape as npReshape
        gradients = self.readkf('GeoOpt','Gradients_CART')
        unitConv = Units.convert(1.0,'a.u.',eUnit) / Units.convert(1.0,'bohr',lUnit)
        #negative of the values in KF are the gradients
        gradients = npAr([ -v * unitConv for v in gradients ])
        nAt = len(gradients)//3
        return npReshape(gradients,(nAt,3))


    def get_energy_decomposition(self, unit='au'):
        """get_energy(unit='au')
        Return a dictionary with energy decomposition terms, expressed in *unit*.

        The following keys are present in the returned dictionary: ``Electrostatic``, ``Kinetic``, ``Coulomb``, ``XC``. The sum of all the values is equal to the value returned by :meth:`get_energy`.
        """
        ret = {}
        ret['Electrostatic'] = self._get_single_value('Energy', 'Electrostatic Energy', unit)
        ret['Kinetic'] = self._get_single_value('Energy', 'Kinetic Energy', unit)
        ret['Coulomb'] = self._get_single_value('Energy', 'Elstat Interaction', unit)
        ret['XC'] = self._get_single_value('Energy', 'XC Energy', unit)
        return ret

    def get_timings(self):
        """get_timings()

        Return a dictionary with timing statistics of the job execution. Returned dictionary contains keys ``cpu``, ``system`` and ``elapsed``. The values are corresponding timings, expressed in seconds.
        """
        last = self.grep_output(' Total Used : ')[-1].split()
        ret = {}
        ret['elapsed'] = float(last[-1])
        ret['system'] = float(last[-3])
        ret['cpu'] = float(last[-5])
        return ret


    def _atomic_numbers_input_order(self):
        """_atomic_numbers_input_order()
        Return a list of atomic numbers, in the input order.
        """
        n = self.readkf('Geometry', 'nr of atoms')
        tmp = self.readkf('Geometry', 'atomtype').split()
        atomtypes = {i+1 : PT.get_atomic_number(tmp[i]) for i in range(len(tmp))}
        atomtype_idx = self.readkf('Geometry', 'fragment and atomtype index')[-n:]
        atnums = [atomtypes[i] for i in atomtype_idx]
        return self.inputOrder(atnums)


    def inputOrder(self,data):
        """_reorder()
        Reorder any iterable data to match the input atom order. Returns a List!
        """
        mapping = self._int2inp()
        return [ d[mapping[i]-1] for i in range(nAt) ]


    def _int2inp(self):
        """_int2inp()
        Get mapping from the internal atom order to the input atom order.
        """
        aoi = self.readkf('Geometry', 'atom order index')
        n = len(aoi)//2
        return aoi[:n]


class ADFJob(SCMJob):
    _result_type = ADFResults
    _command = 'adf'

    def _serialize_mol(self):
        for i,atom in enumerate(self.molecule):
            smb = self._atom_symbol(atom)
            suffix = ''
            if 'adf' in atom.properties and 'fragment' in atom.properties.adf:
                suffix += 'f={fragment} '
            if 'adf' in atom.properties and 'block' in atom.properties.adf:
                suffix += 'b={block}'

            self.settings.input.atoms['_'+str(i+1)] = ('{:>5}'.format(i+1)) + atom.str(symbol=smb, suffix=suffix, suffix_dict=atom.properties.adf)

    def _remove_mol(self):
        if 'atoms' in self.settings.input:
            del self.settings.input.atoms
