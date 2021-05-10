from os.path import join as opj
import numpy as np

from ...core.basejob  import SingleJob
from ...core.settings import Settings
from ...core.errors import PlamsError
from ...mol.molecule  import Molecule
from ...core.results import Results
from ...tools.units import Units

__all__ = ['VASPResults', 'VASPJob']


class VASPResults(Results):
    """A class for VASP results."""
    pass


class VASPJob(SingleJob):
    """
    A class representing a single computational job with `VASP <https://www.vasp.at/>`

    * Set 'ignore_molecule' in ``self.settings`` to disable Molecule handling through ASE.
    * Set 'ignore_potcar' in ``self.settings`` to disable automatic `POTCAR` creation.
    * Set the path to the `POTCAR` files in  ``self.settings.input.potcar`` for automatic `POTCAR` creation.
    * If `POTCAR` files not matching the element symbol should be used, give a translation dict in ``self.settings.input.potcardict``.
    * Settings branch ``input.incar`` is parsed into the `INCAR` file, ``input.xxx`` into the corresponding `XXX` file.
    """
    _command = 'vasp_std'
    _filenames = {'inp':'INCAR', 'run':'$JN.run', 'out':'OUTCAR', 'err': '$JN.err', 'log': '$JN.log'}
    _result_type = VASPResults

    def get_input(self):
        """
        Transform all contents of ``input`` branch of ``settings`` into string.
        """

        def _vaspstr(inp):
            #convert to VASP type string
            if isinstance(inp, bool):
                if inp:
                    return ".TRUE."
                else:
                    return ".FALSE."
            else:
                return str(inp)

        def parse(head, sett):
            ret = ''

            for key in sett:
                if isinstance(sett[key], Settings):
                    raise PlamsError("Nested Settings instance found under input.{}.{}, not supported!".format(head, key))
                elif isinstance(value, list):
                    raise PlamsError("List instance found under input.{}.{}, not supported!".format(head, key))
                else:
                    ret += "{} = {}".format(key,_vaspstr(sett[key]))

            return ret

        use_molecule = ('ignore_molecule' not in self.settings) or (self.settings.ignore_molecule == False)
        use_potcar = ('ignore_potcar' not in self.settings) or (self.settings.ignore_potcar == False)
        if use_molecule:
            self._parsemol()
        if use_potcar:
            self._parsepotcar()

        #we need all first-level keys in uppercase
        tmp = Settings()
        for key in self.settings.input:
            tmp[key.upper()] = self.settings.input[key]

        for item in tmp:
            #POTCAR creation handled above
            if 'POTCAR' in item:
                continue
            #INCAR
            if item == 'INCAR':
                inp = parse(item, tmp[item])
            else:
                with open(opj(self.path,item),'w') as f:
                    f.write(parse(item, tmp[item]))

        return inp

    def _parsemol(self):
        if 'ase' in Molecule._writeformat:
            #ASE has a write function for VASP coordinate files, use that if possible
            filename = opj(self.path, 'POSCAR')
            self.molecule.writease(filename)
        else:
            raise PlamsError('VASP Interface has no builtin Molecule support, install ASE. See Doc for details.')

    def _parsepotcar(self):
        if 'potcar' in self.input:
            elements = [ self.molecule.atoms[0].symbol ]
            for atom in self.molecule.atoms[1:]:
                if not atom.symbol == elements[-1]:
                    elements.append(atom.symbol)
            if 'potcardict' in self.input:
                translate = self.input.potcardict
            else:
                translate = { el: el for el in set(elements)  }
            #open files
            files = [ open(opj(self.input.potcar,translate[el],"POTCAR"),'r') for el in elements ]

            with open(opj(self.path,'POTCAR'), 'w') as f:
                for potcar in files:
                    f.write(potcar.read())
                    potcar.close()
        else:
            raise PlamsError('VASP Interface needs the POTCAR path in self.input.potcar.')


    def get_runscript(self):
        """
        Run VASP

        Overwrite ``self._command`` to change the default VASP Binary.
        """
        ret = self._command
        if self.settings.runscript.stdout_redirect:
            ret += ' >' + self._filename('log')
        ret += '\n\n'
        return ret

    def check(self):
        """
        Look for the normal termination line in output. Note, that does not mean your calculation was successful!
        """
        termination = self.results.grep_output('General timing and')
        return bool(termination)
