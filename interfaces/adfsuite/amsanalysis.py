from ...core.errors import PlamsError
from .scmjob import SCMJob, SCMResults

__all__ = ['AMSAnalysisJob', 'AMSAnalysisResults']

class AMSAnalysisXYResults:
    def __init__(self):
        self.x = []
        self.y = []
        self.x_units = None
        self.y_units = None
        self.y_sigma = [] # stadard deviation for y_values

        self.properties = None
        self.legend = None

    def read_properties (self, kf, sec) :
        """
        Read properties from the KF file
        """
        counter = 0
        properties = {}
        while(1) :
            counter += 1
            try :
                propname = kf.read(sec, 'Property(%i)'%(counter))
            except :
                break
            properties[propname] = kf.read(sec, propname)
        self.properties = properties
        if 'Legend' in properties :
            self.legend = properties['Legend']

class AMSAnalysisResults(SCMResults):
    _kfext = '.kf'
    _rename_map = {'analysis.kf':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('AMSAnalysisResults does not support the get_molecule() method.')

    def get_sections(self) :
        """
        Read the sections available to make xy plots
        """
        if not self._kfpresent():
            raise FileError('File {} not present in {}'.format(self.job.name+self.__class__._kfext, self.job.path))
        sections = self._kf.reader._sections.keys()
        return sections

    def get_xy(self, section='', i=1):
        xy = AMSAnalysisXYResults()

        task = self.job.settings.input.Task
        if section == '' :
            section = task

        # Find the correct section in the KF file
        sections = self.get_sections()
        matches = [s for s in sections if s.lower()==section.lower()+'('+str(i)+')']
        if len(matches) == 0 :
                print ('Sections: ',list(sections))
                raise PlamsError('AMSAnalysisResults.get_xy(section,i): section must be one of the above. You specified "{}"'.format(section))
        sec = matches[0] 

        # Get the data
        x_name = self.readkf(sec, 'x(1)-axis')
        xy.x = self.readkf(sec, x_name) 
        xy.x_units = self.readkf(sec, '%s(units)'%(x_name))

        y_name = self.readkf(sec, 'y-axis')
        xy.y = self.readkf(sec, y_name) 
        xy.y_units = self.readkf(sec, '%s(units)'%(y_name))

        xy.y_sigma = self.readkf(sec, 'sigma')

        xy.read_properties(self._kf, sec)

        return xy

    def get_D(self, i=1):
        """ returns a 2-tuple (D, D_units) from the AutoCorrelation(i) section on the .kf file. """

        sec = 'AutoCorrelation('+str(i)+')'
        D = self.readkf(sec, 'DiffusionCoefficient')
        D_units = self.readkf(sec, 'DiffusionCoefficient(units)')
        return D, D_units


class AMSAnalysisJob(SCMJob):
    """A class for analyzing molecular dynamics trajectories using the ``analysis`` program.


    """
    _result_type = AMSAnalysisResults
    _command = 'analysis'

    def __init__(self, **kwargs):
        SCMJob.__init__(self, **kwargs)

    def _serialize_mol(self):
        pass

    def _remove_mol(self):
        pass

    def check(self):
        try:
            grep = self.results.grep_file('$JN.err', 'NORMAL TERMINATION')
        except:
            return False
        return len(grep) > 0
