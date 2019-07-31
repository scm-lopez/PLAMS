import os
import subprocess

import numpy as np

from scm.plams.tools.units import Units
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)

__all__ = ['CRSResults', 'CRSJob']


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}

    def get_energy(self, unit: str = 'kcal/mol') -> float:
        """Returns the solute solvation energy from an Activity Coefficients calculation."""
        E = self.readkf('ACTIVITYCOEF', 'deltag')[0]
        return Units.convert(E, 'kcal/mol', unit)

    def get_activity_coefficient(self) -> float:
        """Return the solute activity coefficient from an Activity Coefficients calculation."""
        return self.readkf('ACTIVITYCOEF', 'gamma')[0]

    def get_sigma_profile(self, as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, returning a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"σ (e/A**2)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self._get_sigma('SIGMAPROFILE', as_df=as_df)
        except KeyError:
            return self._get_sigma('PURESIGMAPROFILE', as_df=as_df)

    def get_sigma_potential(self, unit: str = 'kcal/mol', as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, expressed in *unit* and return a dictionary of Numpy Arrays.

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
            return self._get_sigma('SIGMAPOTENTIAL', unit=unit, as_df=as_df)
        except KeyError:
            return self._get_sigma('PURESIGMAPOTENTIAL', unit=unit, as_df=as_df)

    def _get_sigma(self, section: str, unit: str = 'kcal/mol', as_df: bool = False) -> dict:
        """Grab all values of sigma and the sigmapotential/profile;
        combine them into a dictionary or pandas dataframe.

        """
        subsection = 'mu' if 'POTENTIAL' in section else 'profil'
        ret = self._sigma_y(section, subsection, unit)
        ret['σ (e/A**2)'] = np.array(self.readkf(section, 'chdval'))

        # Return the dictionary of sigma profiles/potentials
        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section)

    def _sigma_y(self, section: str, subsection: str, unit: str) -> dict:
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

    @staticmethod
    def _dict_to_df(sigm_dict: dict, section: str) -> 'pandas.DataFrame':
        """Attempt to convert a dictionary, produced by :meth:`._get_sigma` into a DataFrame."""
        try:
            import pandas as pd
            index = pd.Index(sigm_dict.pop('σ (e/A**2)'), name='σ (e/A**2)')
            df = pd.DataFrame(sigm_dict, index=index)
            df.columns.name = section.lower()
            return df
        except ImportError as ex:
            raise ex.__class__("CRSResults._get_sigma: as_df=True requires the 'pandas' package")


class CRSJob(SCMJob):
    """A |SCMJob| subclass intended for running COSMO-RS jobs."""
    _command = 'crs'
    _result_type = CRSResults

    def __init__(self, **kwargs) -> None:
        """Initialize a :class:`CRSJob` instance."""
        super().__init__(**kwargs)
        self.settings.ignore_molecule = True

    @staticmethod
    def cos_to_coskf(filename: str) -> str:
        """Convert a .cos file into a .coskf file with the :code:`$ADFBIN/cosmo2kf` command.
        Returns the filename of the new .coskf file.

        """
        filename_out = filename + 'kf'
        adfbin = os.environ['ADFBIN']

        args = [os.path.join(adfbin, 'cosmo2kf'), filename, filename_out]
        subprocess.run(args)
        return filename_out
