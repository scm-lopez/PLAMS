from ...core.errors import PlamsError
from .scmjob import SCMJob, SCMResults

__all__ = ['AMSAnalysisJob', 'AMSAnalysisResults','convert_to_unicode']

class AMSAnalysisXYResults:
    def __init__(self):
        self.x = []
        self.y = []
        self.x_units = None
        self.y_units = None
        self.x_name = None
        self.y_name = None
        self.y_sigma = [] # stadard deviation for y_values

        self.properties = None
        self.legend = None

    def read_data (self, kf, sec) :
        """
        Read the xy data for a section from the kf file
        """
        xkey = 'x(1)-axis'
        x_name = kf.read(sec, '%s(label)'%(xkey))
        self.x = kf.read(sec, xkey)
        self.x_name = convert_to_unicode(x_name)
        self.x_units = convert_to_unicode(kf.read(sec, '%s(units)'%(xkey)))

        ykey = 'y-axis'
        y_name = kf.read(sec, '%s(label)'%(ykey))
        self.y = kf.read(sec, ykey)
        self.y_name = convert_to_unicode(y_name)
        self.y_units = convert_to_unicode(kf.read(sec, '%s(units)'%(ykey)))

        self.y_sigma = kf.read(sec, 'sigma')

        self.read_properties(kf, sec)

    def read_properties (self, kf, sec) :
        """
        Read properties from the KF file
        """
        counter = 0
        properties = {}
        while(1) :
            counter += 1
            try :
                propname = kf.read(sec, 'Property(%i)'%(counter)).strip()
            except :
                break
            properties[propname] = kf.read(sec, propname)
            if isinstance(properties[propname],str) :
                properties[propname] = properties[propname].strip()

        # Now set the instance variables
        self.properties = properties
        if 'Legend' in properties :
            self.legend = properties['Legend']

    def write_plot (self, outfilename=None) :
        """
        Print this plot to a text file
        """
        parts = []
        for propname,prop in self.properties.items() :
            parts.append('%-30s %s\n'%(propname, prop))
        x_name = '%s(%s)'%(self.x_name,self.x_units)
        y_name = '%s(%s)'%(self.y_name,self.y_units)
        parts.append('%30s %30s %30s\n'%(x_name,y_name,'sigma'))
        for x,y,s in zip(self.x,self.y,self.y_sigma) :
            parts.append('%30.10e %30.10e %30.10e\n'%(x,y,s))
        block = ''.join(parts)

        if outfilename is not None :
            outfile = open(outfilename,'w',encoding='utf8')
            outfile.write(block)
            outfile.close()

        return block
        

class AMSAnalysisResults(SCMResults):
    _kfext = '.kf'
    _rename_map = {'plot.kf':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('AMSAnalysisResults does not support the get_molecule() method.')

    def get_sections(self) :
        """
        Read the sections available to make xy plots
        """
        if not self._kfpresent():
            raise FileError('File {} not present in {}'.format(self.job.name+self.__class__._kfext, self.job.path))
        if self._kf.reader._sections is None :
            self._kf.reader._create_index()
        sections = self._kf.reader._sections.keys()
        return sections

    def get_xy(self, section='', i=1):
        xy = AMSAnalysisXYResults()

        task = self.job.settings.input.Task
        if section == '' :
            section = task

        # Find the correct section in the KF file
        sections = self.get_sections()
        matches = [s for s in sections if s.lower()==section.lower()+'(%i)'%(i)]
        if len(matches) == 0 :
                print ('Sections: ',list(sections))
                raise PlamsError('AMSAnalysisResults.get_xy(section,i): section must be one of the above. You specified "{}"'.format(section))
        sec = matches[0] 

        # Get the data
        xy.read_data(self._kf,sec)

        return xy

    def write_all_plots (self) :
        """
        Get a list of all the plots created by the analysis job
        """
        sections = self.get_sections()
        for section in sections :
            name_part = section.split('(')[0]
            num_part = int(section.split('(')[1].split(')')[0])
            outfilename = '%s_%i.dat'%(name_part,num_part)
            xy = self.get_xy(name_part,num_part)
            xy.write_plot('%s'%(outfilename))

    def get_D(self, i=1):
        """ returns a 2-tuple (D, D_units) from the AutoCorrelation(i) section on the .kf file. """

        # If there are multiple, it will read the first one
        sections = [sec for sec in self.get_sections() if 'Integral' in sec]
        if len(sections) < i : 
            return None,None
        section = sections[i-1]
        plot = self.get_xy(section.split('(')[0],i)
        if not 'DiffusionCoefficient' in plot.properties.keys() :
            return None, None

        D = plot.properties['DiffusionCoefficient']
        D_units = plot.y_units
        return D, D_units


class AMSAnalysisJob(SCMJob):
    """A class for analyzing molecular dynamics trajectories using the ``analysis`` program.


    """
    _result_type = AMSAnalysisResults
    _command = 'analysis'
    _subblock_end = 'end'

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

def convert_to_unicode (k) :
    """
    Convert a string with ascii symbols representing unicode symbols

    Example k: 'abc\\u03c9def'
    """
    parts = k.split('\\u')
    # Collect the hexadecimals
    symbols = [chr(int(part[:4],16)) for part in parts[1:]]
    # Now repair the parts
    parts = parts[:1] + [''.join([s,part[4:]]) for s,part in zip(symbols,parts[1:])]
    key = ''.join(parts)

    return key
