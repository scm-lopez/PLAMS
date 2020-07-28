"""
Run RASPA2 with PLAMS
Contributed by Patrick Melix

See the Documentation for an Example
"""
from os.path import join as opj

from ...core.basejob import SingleJob
from ...core.settings import Settings
from ...core.results import Results
from ...tools.units import Units

__all__ = ['RaspaJob', 'RaspaResults']


class RaspaResults(Results):
    """A Class for handling RASPA Results.

    Note that the Output file is guessed in ``self.collect``. Do not rely on it to point to the
    file you need and but check ``self._filenames['out']``!
    """

    def collect(self):
        """Try to Guess the Main Output File.

        Set it to whatever comes first in the 'Output/System_0' folder.
        """
        Results.collect(self)
        if 'out' not in self.job._filenames:
            files = [f for f in self.files if 'Output/System_0' in f]
            if files:
                self.job._filenames['out'] = files[0]






class RaspaJob(SingleJob):
    """A single job with RASPA2 class.

    Requires the ASE Python package.

    Files that should be copied to the job directory can be passed as a dictionary using the ``copy``
    argument. The dict should have the following layout: *dict[<filename in jobdir>] = <path to file>*.

    **Molecule parsing is not yet supported!**
    """
    _result_type = RaspaResults
    _filenames = {'inp':'simulation.input', 'run':'$JN.run', 'err': '$JN.err', 'stdout': '$JN.out'}


    def __init__(self, copy=None, **kwargs):
        SingleJob.__init__(self, **kwargs)
        self.copy_files = copy


    def _get_ready(self):
        """Copy files to execution dir if self.copy_files is set."""
        SingleJob._get_ready(self)
        if self.copy_files:
            for filename, path in self.copy_files.items():
                dest = opj(self.path, filename)
                shutil.copy(path, dest)
        return


    def get_input(self):
        """Transform all contents of ``settings.input`` branch into string with blocks, keys and values.

        RASPA2 does not interpret indentation but rather certain keys automatically start a subsection.
        Subsections are never closed but *automatically* ended if the next key-value pair does not fit into it.
        So here we only need to care about readability of the input file. The logic is automatically preserved
        by using the nested PLAMS ``settings`` logic.
        """

        def parse(key, value, indent=''):
            if value is True:
                value = 'yes'
            elif value is False:
                value = 'no'

            ret = indent + key

            if isinstance(value, Settings):
                #reserved keywords are not part of the input file
                if '_h' in value:
                    ret += ' ' + value['_h']
                ret += '\n'

                for el in value:
                    if el != '_h':
                        ret += parse(el, value[el], indent+'  ')
            else:
                ret += ' ' + str(value)
            ret += '\n'
            return ret

        inp = ''
        for item in self.settings.input:
            inp += parse(item, self.settings.input[item])

        return inp


    def get_runscript(self):
        """$RASPA_DIR has to be set in the environment!"""
        ret = 'export DYLD_LIBRARY_PATH=${RASPA_DIR}/lib'
        ret += '\nexport LD_LIBRARY_PATH=${RASPA_DIR}/lib'
        ret += '\n$RASPA_DIR/bin/simulate '
        if self.settings.runscript.stdout_redirect:
            ret += ' >' + self._filename('stdout')
        ret += '\n\n'
        return ret


    def check(self):
        """Look for the normal termination signal in output."""
        s = self.results.grep_output('Simulation finished on')
        return len(s) == 0
