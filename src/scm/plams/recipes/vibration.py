from collections import OrderedDict
from ..core.results import Results
from ..core.basejob import MultiJob
from ..core.functions import log
from ..interfaces.adfsuite.adf import ADFJob
from numpy import array as npa
from pickle import dump as pickleDump
from os.path import join as osPJ
from os.path import basename as osPB

__all__ = ['VibrationsResults' ,'VibrationsJob']

try:
    from ase.vibrations import Vibrations as aseVib
    from ..tools.ase import *
except ImportError:
    __all__ = []


class VibrationsResults(Results):
    """
    A Class for handling numerical Vibrational analysis Results.

    Use ``get_ASEVib()`` to access the ``ase.vibrations.Vibrations`` object. Using e.g. ``self.get_ASEVib().summary()`` the results of the ASE frequency calculation are accessible.
    """
    def get_ASEVib(self):
        """Returns a ``ase.vibrations.Vibrations`` object"""
        return self._vib


class VibrationsJob(MultiJob):
    """
    Class for calculating numerical Frequencies in parallel using arbitrary interfaces.
    This is achieved using the Vibrations class from ASE.

    *   ``molecule`` -- |Molecule| object (most propably in the chosen Methods optimized state)
    *   ``settings`` -- |Settings| instance for all Single-Point jobs to be run. Don't forget reference to a restart file if you want to save a lot of computer time!
    *   ``jobType`` -- |Job| Class you want to use.
    *   ``get_gradients`` -- Function name to retrieve gradients of the |Results| object of your chosen ``jobType.results``. Must take options ``eUnit='eV'`` and ``lUnit='Angstrom'``.
    *   ``reorder`` -- Function name of ``jobType.results`` to reorder gradients to input order, set to ``None`` if not applicable.
    *   ``aseVibOpt`` -- Options for ``ase.vibrations.Vibrations.__init__``

    The ``self.__init__()`` method calls on ase.vibrations.Vibrations.run to create all displacements, then adds them as children.
    After running, the ``self.postrun()`` method reads the gradients from all single-points and saves them according to the ASE scheme.
    Finally the ASE Vibrations object is made accessible in ``self.results.get_ASEVib`` for the retrieval of the ASE results.

    """
    _result_type = VibrationsResults

    def __init__(self, molecule, settings, jobType=ADFJob, get_gradients='get_gradients', reorder='inputOrder', aseVibOpt={}):
        MultiJob.__init__(self)
        self.molecule = molecule
        self.settings = settings
        self.jobType = jobType
        self.get_gradient = get_gradients
        self.reorder = reorder
        self.aseVibOpt = aseVibOpt

        self.get_grad = getattr(self.jobType._result_type, get_gradients)
        if reorder != None:
            self.reorder = getattr(self.jobType._result_type, reorder)
        else:
            self.reorder = lambda x: x


    def new_children(self):
        #Create displaced molecules with ASE
        #Include our workingdir in the name, this results in the ASE .pckl files to be saved there
        #kinda hacky but works

        if len(self.children) == 0:
            add = []
            path = self.path
            self._vib = aseVib(toASE(self.molecule), name=osPJ(path,'plams.vib'), **self.aseVibOpt)
            aseMols = self._vib.run(export='return')
            for key, item in aseMols.items():
                add.append(self.jobType(molecule=fromASE(item), name=osPB(key).replace('.pckl',''), settings=self.settings))
            return add

        else:
            return None


    def postrun(self):
        #get gradients and save them to the pickle files for ASE
        path = self.path
        for child in self.children:
            res = child.results
            name = child.name
            # don't rely on gradients beeing numpy arrays
            f = [ npa(vec) for vec in self.get_grad(res, eUnit='eV', lUnit='Angstrom') ]
            force = -1.0 * npa(self.reorder(res, f))
            filename = osPJ(path, name+'.pckl')
            with open(filename, 'wb') as f:
                pickleDump(force, f, protocol=2)

        self.results._vib = self._vib
