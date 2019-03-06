from ..interfaces.adfsuite.adf import ADFJob
from ..core.functions import log

__all__ = ['ADFNBOJob']

class ADFNBOJob(ADFJob):

    def prerun(self):
        s = self.settings.input
        s[s.find_case('fullfock')] = True
        s[s.find_case('aomat2file')] = True
        s[s.find_case('symmetry')] = 'NoSym'
        basis = s.find_case('basis')
        core = s[basis].find_case('core')
        s[basis][core] = 'None'
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
        self.settings.runscript.post = '$ADFBIN/adfnbo <<eor\n' + '\n'.join(adfnbo_input) + '\neor\n\n$ADFBIN/gennbo6 FILE47\n'
