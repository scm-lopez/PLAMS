import os
import subprocess

import numpy as np

from scm.plams.core.basejob import SingleJob, Job, Results, ResultsError
from scm.plams.tools.units import Units
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)

__all__ = ['CRSResults', 'CRSJob']


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}

    def get_energy(self, unit='kcal/mol') -> float:
        """Returns the solute solvation energy from an Activity Coefficients calculation."""
        E = self.readkf('ACTIVITYCOEF', 'deltag')[0]
        return Units.convert(E, 'kcal/mol', unit)

    def get_activity_coefficient(self) -> float:
        """Return the solute activity coefficient from an Activity Coefficients calculation."""
        return self.readkf('ACTIVITYCOEF', 'gamma')[0]

    def get_sigma_profile(self, as_df=False) -> dict:
        r"""Grab all sigma profiles, returning a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"'σ (e/A**2)'"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self.get_sigma('SIGMAPROFILE', as_df=as_df)
        except KeyError:
            return self.get_sigma('PURESIGMAPROFILE', as_df=as_df)

    def get_sigma_potential(self, unit='kcal/mol', as_df=False) -> dict:
        r"""Grab all sigma profiles, expressed in *unit*, returning a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"'σ (e/A**2)'"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self.get_sigma('SIGMAPOTENTIAL', unit=unit, as_df=as_df)
        except KeyError:
            return self.get_sigma('PURESIGMAPOTENTIAL', unit=unit, as_df=as_df)

    def get_sigma(self, section, unit='kcal/mol', as_df=False) -> dict:
        """Grab all values of sigma and the sigmapotential/profile;
        combine them into a dictionary or pandas dataframe.

        """

        def _sigma_y(section, subsection, unit) -> dict:
            """Construct all values for the y-axis"""
            # Use filenames as keys
            _filenames = self.readkf(section, 'filename').split()
            filenames = [_filenames] if not isinstance(_filenames, list) else _filenames
            keys = [os.path.split(key)[-1] for key in filenames]

            # Use sigma profiles/potentials as values
            values = np.array(self.readkf(section, subsection))
            values *= Units.conversion_ratio('kcal/mol', unit)
            values.shape = len(keys), len(values) // len(keys)

            ret = dict(zip(keys, values))
            if 'PURE' not in section:  # Add a final key for the solvent mixture
                ret['Total'] = np.array(self.readkf(section, subsection + 'tot'))
            return ret

        subsection = 'mu' if 'POTENTIAL' in section else 'profil'
        ret = _sigma_y(section, subsection, unit)
        ret['σ (e/A**2)'] = np.array(self.readkf(section, 'chdval'))

        # Return the dictionary of sigma profiles/potentials
        if not as_df:
            return ret

        # Attempt to convert the to-be returned dictionary into a Pandas DataFrame
        try:
            import pandas as pd
            index = pd.Index(ret.pop('σ (e/A**2)'), name='σ (e/A**2)')
            df = pd.DataFrame(ret, index=index)
            df.columns.name = section.lower()
            return df
        except ModuleNotFoundError:
            raise ModuleNotFoundError("get_sigma: as_df=True requires the 'Pandas' package")


class CRSJob(SCMJob):
    """A |SCMJob| subclass intended for running COSMO-RS jobs."""
    _command = 'crs'
    _result_type = CRSResults

    def __init__(self, **kwargs) -> None:
        """Initialize a :class:`CRSJob` instance."""
        SingleJob.__init__(self, **kwargs)
        self.settings.ignore_molecule = True

    def _unpack_results(self, results) -> str:
        """Search for .cos, .coskf, .crskf or .t21 files in *results*.

        Returns the path+filename of the first matching file.
        Raises a :exc:`ResultsError` if no such files can be found in *results*.

        """
        extensions = ('.crskf', '.coskf', '.cos', '.t21')

        # Search for .crskf, .coskf or .t21 files
        for file in results.files:
            for ext in extensions:
                if ext not in file:
                    continue

                # A match has been found!
                path = results.job.path
                ret = os.path.join(path, file)
                if ext == '.cos':
                    self._cos_to_coskf(ret)
                    ret += 'kf'  # Convert .cos into .coskf
                return ret

        # No .crskf, .coskf, .cos or .t21 file found, raise a ResultsError
        raise ResultsError('_get_crs_file: A Results instance was passed to CRSJob.settings'
                           ' which does not contain any .cos, .crskf, .coskf or .t21 files')

    @staticmethod
    def _cos_to_coskf(filename) -> str:
        """Convert a .cos file into a .coskf file with the :code:`$ADFBIN/cosmo2kf` command."""
        output = filename + 'kf'
        adfbin = os.environ['ADFBIN']
        args = [os.path.join(adfbin, 'cosmo2kf'), filename, output]
        subprocess.run(args)

    def run(self, jobrunner=None, jobmanager=None, **kwargs) -> CRSResults:
        """Run the job using *jobmanager* and *jobrunner* (or defaults, if ``None``).

        Functions exactly the same as the default :meth:`.Job.run` method barring one exception:
        If a :class:`.Results` instance is assigned to the ``['input']['compound']`` key in
        :attr:`CRSJob.settings` then this method will attempt to extract the path+filename of
        a .cos, .coskf, .crskf or .t21 file from aforementiond :class:`.Results` instance.
        A :exc:`ResultsError` is raised if a :class:`.Results` instance is found that does
        not contain any .cos, .coskf, .crskf or .t21 files.

        """
        # Check if the 'compound' key is present in the settings
        s = self.settings.input
        key = s.find_case('compound')
        if key not in s:
            return Job.run(self, jobrunner=None, jobmanager=None, **kwargs)
        value = s[key]

        # Value is a list of settings; search for any Results instances and convert to a string
        if isinstance(value, list):
            for item in value:
                if not isinstance(item._h, Results):
                    continue
                item._h = self._unpack_results(item._h)

        # Value is Results instance; convert it into a string
        elif isinstance(value._h, Results):
            value._h = self._unpack_results(value._h)

        return Job.run(self, jobrunner=None, jobmanager=None, **kwargs)
