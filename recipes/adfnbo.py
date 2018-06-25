from ..interfaces.adfsuite.adf import ADFJob
from ..core.functions import log

__all__ = ['ADFNBOJob']

class ADFNBOJob(ADFJob):

    def prerun(self):
        s = self.settings.input
        s[s.find_case('fullfock')] = True
        s[s.find_case('aomat2file')] = True
        s[s.find_case('symmetry')] = 'NOSYM'
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

        self.settings.runscript.post =  """
$ADFBIN/adfnbo <<eor
write
fock
spherical
end input
eor

$ADFBIN/gennbo6 FILE47
"""
