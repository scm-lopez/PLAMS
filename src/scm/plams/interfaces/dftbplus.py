"""
Run DFTB+ with plams
v0.1 by Patrick Melix

add 'from .dftbplusjob import *' to $ADFHome/scripting/plams/src/scm/plams/___init___.py to make this work
add this file to $ADFHome/scripting/plams/src/scm/plams/

see dftb+_example.plms for an example
"""
from __future__ import unicode_literals


from .basejob import SingleJob
from .settings import Settings
from .results import Results
from .utils import Units
from .basemol import Molecule
from .errors import PlamsError
import os

__all__ = ['DFTBPlusJob']

class DFTBPlusResults(Results):
    """A Class for DFTB+ Results"""
    _outfile = 'detailed.out'
    _xyzout = 'geo_end.xyz'

    def get_molecule(self):
        #read molecule coordinates
        try:
            mol = Molecule(filename=self[self._xyzout])
        except:
            mol = Molecule()
        return mol

    def get_energy(self, string='Total energy',unit='au'):
        #get the energy given in the output with description "<string>:"
        try:
            energy = float(self.grep_file(self._outfile, pattern=string+':')[0].split()[2])
            energy = Units.convert(energy, 'au', unit)
        except:
            energy = float('nan')
        return energy

    def get_atomic_charges(self):
        #returns dictonary with atom numbers and their charges
        try:
            atomic_charges = {}
            string = self.awk_file(self._outfile,script='/Net atomic charges/{do_print=1} NF==0 {do_print=0 }do_print==1 {print}')
            for line in string:
                if string[0] == line or string[1] == line:
                    continue
                l = line.split()
                atomic_charges[l[0]] = float(l[1])
        except:
            atomic_charges = {}
        return atomic_charges



class DFTBPlusJob(SingleJob):
    """A class representing a single computational job with DFTB+.
       Only supports molecular coordinates, no support for lattice yet"""
    _result_type = DFTBPlusResults
    _filenames = {'inp':'dftb_in.hsd', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}


    def get_input(self):
      """Transform all contents of ``setting.input`` branch into string with blocks, keys and values.

      Automatic handling of ``molecule`` can be disabled with ``settings.ignore_molecule = True``.
      """
      def parse(key, value, indent=''):
          if value is True:
              value = ''

          ret = indent + key
          if key != '' and value != '':
              ret += ' ='

          if isinstance(value, Settings):
              if '_h' in value:
                  ret += ' ' + value['_h']
              ret += ' {\n'

              i = 1
              while ('_'+str(i)) in value:
                  ret += parse('', value['_'+str(i)], indent+'  ')
                  i += 1

              for el in value:
                  if not el.startswith('_'):
                      ret += parse(el, value[el], indent+'  ')
              ret += indent + '}'
          else:
              ret += ' ' + str(value)
          ret += '\n'
          return ret

      inp = ''
      use_molecule = ('ignore_molecule' not in self.settings) or (self.settings.ignore_molecule == False)
      if use_molecule:
          self._parsemol()

      for item in self.settings.input:
          inp += parse(item, self.settings.input[item])

      if use_molecule:
          self._removemol()
      return inp


    def _parsemol(self):
        atom_types = {}
        n = 1
        atoms_line = ''
        for atom in self.molecule:
            if atom.symbol not in atom_types:
                atoms_line += atom.symbol + ' '
                atom_types[atom.symbol] = n
                n += 1

        self.settings.input.geometry._h = 'GenFormat'
        self.settings.input.geometry._1 = '%i C'%len(self.molecule)
        self.settings.input.geometry._2 = atoms_line
        self.settings.input.geometry._3 = ''

        for i,atom in enumerate(self.molecule):
            self.settings.input.geometry['_'+str(i+4)] = ('%5i'%(i+1)) + atom.str(symbol=str(atom_types[atom.symbol]))


    def _removemol(self):
        if 'geometry' in self.settings.input:
            del self.settings.input.geometry


    def get_runscript(self):
        #dftb+ has to be the bin to run
        ret = 'dftb+ '
        if self.settings.runscript.stdout_redirect:
            ret += ' >' + self._filenames('out')
        ret += '\n\n'
        return ret


    def check(self):
        s = self.results.grep_output('ERROR!')
        return len(s) == 0

