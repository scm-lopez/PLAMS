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
    """A Class for handling numerical Vibrational analysis Results"""
    pass

class VibrationsJob(MultiJob):
    """
    Class for calculating numerical Frequencies in parallel using arbitrary interfaces.
    This is achieved using the Vibrations class from ASE.

    *   ``molecule`` -- Molecule (most propably in the chosen Methods optimized state)
    *   ``settings`` -- |Settings| instance for all Single-Point jobs to be run. Don't forget reference to a restart file if you want to save a lot of computer time!
    *   ``jobType`` -- |Job| Class you want to use
    *   ``get_gradients`` -- Function name to retrieve gradients of the |Result| object of your chosen jobType.results. Must take options eUnit='eV' and lUnit='Angstrom'
    *   ``reorder`` -- Function name of jobType.results to reorder gradients to input order, set to ``None`` if not applicable.
    *   ``aseVibOpt`` -- Options for ase.vibrations.Vibrations.__init__

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

    def prerun(self):
        #Create displaced molecules with ASE
        #Include our workingdir in the name, this results in the ASE .pckl files to be saved there
        #kinda hacky but works
        path = self.path
        self.vib = aseVib(toASE(self.molecule), name=osPJ(path,'plams.vib'), **self.aseVibOpt)
        aseMols = self.vib.run(export='return')
        for key, item in aseMols.items():
            self.children.append(self.jobType(molecule=fromASE(item), name=osPB(key).replace('.pckl',''), settings=self.settings))


    def postrun(self):
        #get gradients and save them to the pickle files for ASE
        path = self.path
        for child in self.children:
            res = child.results
            name = child.name
            force = -1.0 * npa(self.reorder(res, self.get_grad(res, eUnit='eV', lUnit='Angstrom')))
            filename = osPJ(path, name+'.pckl')
            with open(filename, 'wb') as f:
                pickleDump(force, f, protocol=2)

        self.results.vib = self.vib
