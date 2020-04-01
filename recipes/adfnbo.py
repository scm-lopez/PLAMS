from ..interfaces.adfsuite.ams import AMSJob
from ..core.functions import log
from ..core.settings import ig
import os

__all__ = ['ADFNBOJob']

class ADFNBOJob(AMSJob):

    def prerun(self):
        s = self.settings.input.ADF
        s[ig('fullfock')] = True
        s[ig('aomat2file')] = True
        s[ig('symmetry')] = 'NoSym'
        s[ig('basis')][ig('core')] = 'None'
        save = s.find_case('save')
        if save in s:
            if isinstance(s.save, str):
                s.save += ' TAPE15'
            elif isinstance(s.save, list):
                s.save.append('TAPE15')
            else:
                log("WARNING: 'SAVE TAPE15' could not be added to the input settings of {}. Make sure (thisjob).settings.input.save is a string or a list.".format(self.name), 1)
        else:
            s[save] = 'TAPE15'

        if isinstance(self.settings.adfnbo, list):
            adfnbo_input = self.settings.adfnbo
        else:
            adfnbo_input = ['write', 'spherical', 'fock']
            log('WARNING: (thisjob).settings.adfnbo should be a list. Using default settings: write, fock, spherical', 1)

        self.settings.runscript.post = 'cp "'+os.path.join(self.path,'adf.rkf')+'" TAPE21\n' '$ADFBIN/adfnbo <<eor\n' + '\n'.join(adfnbo_input) + '\neor\n\n$ADFBIN/gennbo6 FILE47\n'
