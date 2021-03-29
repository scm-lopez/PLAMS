#!/usr/bin/env python

from ..core.settings import Settings
from ..core.errors import PlamsError
from ..tools.kftools import KFFile
from ..interfaces.adfsuite.amsanalysis import AMSAnalysisJob
from .rkffile import RKFTrajectoryFile
from .rkfhistoryfile import RKFHistoryFile

__all__ = ['Trajectory']

class Trajectory :
        """
        Class representing an AMS trajectory

        It creates molecule objects along the trajectory, and it performs analysis
        """
        def __init__ (self, filenames) :
                """
                Initiates instance of an AMS trajectory object

                * ``filenames`` -- List of filepaths of RKF trajectory files or single filepath for RKF file

                Note: The corresponding file needs to remain on disc
                """
                if isinstance(filenames,str):
                        filenames = [filenames]

                # Ceck is RKFHistoryFile needs to be used
                classes = []
                for fn in filenames :
                        if not fn.split('.')[-1] == 'rkf' :
                                raise PlamsError('Files need to be RKF')
                        kf = KFFile(fn)
                        if kf.reader is None :
                                raise PlamsError('RKF file %s not found'%(fn))
                        kf.reader._create_index()
                        if not 'History' in kf.reader._sections.keys() :
                                raise PlamsError('%s is not a trajectory file'%(fn))
                        if 'SystemVersionHistory' in kf.reader._sections.keys() :
                                classes.append(RKFHistoryFile)
                        else :
                                classes.append(RKFTrajectoryFile)

                # Store the file objects and associated data for each file
                self.files = [classes[i](fn) for i,fn in enumerate(filenames)]
                self.molecules = [rkf.get_plamsmol() for rkf in self.files]
                self.lengths = [len(rkf) for rkf in self.files]

        def __len__ (self) :
                """
                Magic method that returns the length of the trajectory
                """
                return sum(self.lengths)

        def __iter__ (self) :
                """
                Iterates over molecule objects
                """
                for irkf,rkf in enumerate(self.files) :
                        mol = self.molecules[irkf]
                        for i in range(self.lengths[irkf]) :
                                print (i) 
                                crd,cell = rkf.read_frame(i,molecule=mol)
                                yield mol.copy()

        def __getitem__ (self, i) :
                """
                Returns Molecule object
                """
                irkf, istep = self._get_filenum_and_stepnum(i)
                mol = self.molecules[irkf].copy()
                crd,cell = self.files[irkf].read_frame(istep,molecule=mol)
                self.files[irkf].rewind()
                return mol

        def run_analysis (self, settings, range=None) :
                """
                Calls the AMS analysis tool behind the scene

                * ``settings`` -- PLAMS Settings object
                                  Example :
                                  settings = Settings()
                                  settings.input.Task = 'AutoCorrelation' 
                                  settings.input.AutoCorrelation.Property = 'Velocities'
                                  settings.input.AutoCorrelation.MaxStep = 2000
                * ``range``    -- Not implemented yet
                """
                trajecsettings = []
                for ikf,rkf in enumerate(self.files) :
                        s = Settings()
                        s.Trajectory.KFFilename = self.files[0].file_object.path
                        #s.Trajectory.Range = '1 10000'
                        trajecsettings.append(s)
                settings.input.TrajectoryInfo = trajecsettings

                job = AMSAnalysisJob(settings=settings)
                result = job.run()
                plots = result.get_all_plots()
                return plots

        #################
        # Private methods
        #################

        def _get_filenum_and_stepnum (self, i) :
                """
                Connects a filenumber and a stepnumber to index i
                """
                irkf,istep = (-1,-1)
                counter = 0
                for il,length in enumerate(self.lengths) :
                        if i < counter+length :
                                irkf = il
                                istep = i-counter
                                break
                        counter += length
                if irkf < 0 :
                        raise KeyError('Index out of range')
                return irkf, istep
