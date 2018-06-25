import subprocess

from ..core.basejob  import SingleJob
from ..core.settings import Settings



class Cp2kJob(SingleJob):
    """
    A class representing a single computational job with `CP2K <https://www.cp2k.org/>`
    """
    def get_input(self):
        """
        Transform all contents of ``input`` branch of ``settings`` into string
        with blocks, subblocks, keys and values.
        """

        _reserved_keywords = ["KIND", "XC", "JOB", "AT_SET", "AT_INCLUDE", "AT_IF"]

        def parse(key, value, indent=''):
            ret = ''
            key = key.upper()
            print("Keys: ", key)
            if isinstance(value, Settings):
                if not any(k in key for k in _reserved_keywords):
                    ret += '{}&{}\n'.format(indent, key)
                    for el in value:
                        ret += parse(el, value[el], indent + '  ')
                    ret += '{}&END\n'.format(indent)

                elif key == "XC":
                    ret += '{}&{}\n'.format(indent, key)
                    for el in value:
                        if el.upper() == "XC_FUNCTIONAL":
                            v = value[el]
                            ret += '  {}&XC_FUNCTIONAL {}\n'.format(indent, v)
                            ret += '  {}&END\n'.format(indent)
                        else:
                            ret += parse(el, value[el], indent + '  ')
                    ret += '{}&END\n'.format(indent)

                elif "KIND" in key:
                    for el in value:
                        ret += '{}&{}  {}\n'.format(indent, key, el.upper())
                        for v in value[el]:
                            ret += parse(v, value[el][v], indent + '  ')
                        ret += '{}&END\n'.format(indent)
                elif "JOB" in key:
                    work_dirs = value['directories']
                    job_names = value['input_file_names']
                    job_ids = value['job_ids']

                    for k, (jobID, name, wd) in enumerate(zip(
                            job_ids, job_names, work_dirs)):
                        ret += '{}&JOB\n'.format(indent)
                        if k > 0:
                            ret += '  {}DEPENDENCIES {}\n'.format(indent, jobID - 1)
                        ret += '  {}DIRECTORY {}\n'.format(indent, wd)
                        ret += '  {}INPUT_FILE_NAME {}\n'.format(indent, name)
                        ret += '  {}JOB_ID {}\n'.format(indent, jobID)
                        ret += '{}&END JOB\n\n'.format(indent)
                elif "AT_SET" in key:
                    var, val = tuple(value.items())[0]
                    ret += '@SET {} {}\n'.format(var, val)

                elif "AT_IF" in key:
                    pred, branch = tuple(value.items())[0]
                    ret += '{}@IF {}\n'.format(indent, pred)
                    for k, v in branch.items():
                        ret += parse(k, v, indent + '  ')
                    ret += '{}@ENDIF\n'.format(indent)

            elif key == "AT_INCLUDE":
                ret += '@include {}\n'.format(value)

            elif isinstance(value, list):
                for el in value:
                    ret += parse(key, el, indent)

            elif value is '' or value is True:
                ret += '{}{}\n'.format(indent, key)
            else:
                ret += '{}{}  {}\n'.format(indent, key, str(value))
            return ret

        inp = ''
        for item in self.settings.input:
            inp += parse(item, self.settings.input[item]) + '\n'

        return inp

    def get_runscript(self):
        """
        Run parallel version of Cp2k using srun.
        """
        # try to cp2k using srun
        try:
            subprocess.run(["srun", "--help"], stdout=subprocess.DEVNULL)
            ret = 'srun cp2k.popt'
        except OSError:
            ret = 'cp2k.popt'

        ret += ' -i {} -o {}'.format(self._filename('inp'), self._filename('out'))

        return ret

    def check(self):
        """
        Look for the normal termination signal in Cp2k output
        """
        s = self.results.grep_output("PROGRAM STOPPED IN")
        return len(s) > 0
