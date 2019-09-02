import os
import subprocess

import numpy as np
try:
    import pandas as pd
    PANDAS = True
except ImportError:
    PANDAS = False

from scm.plams.tools.units import Units
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)

__all__ = ['CRSResults', 'CRSJob']


class CRSResults(SCMResults):
    """A |SCMResults| subclass for accessing results of |CRSJob|."""
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}

    def readarray(self, section: str, subsection: str, **kwargs) -> np.ndarray:
        """Read data from *section*/*subsection* of the main KF file and return as NumPy array.

        All additional provided keyword arguments will be passed onto the numpy.array_ function.

        .. _numpy.array: https://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html
        """
        return np.array(self.readkf(section, subsection), **kwargs)

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
            return self._get_sigma('SIGMAPROFILE', 'profil', as_df=as_df)
        except KeyError:
            return self._get_sigma('PURESIGMAPROFILE', 'profil', as_df=as_df)

    def get_sigma_potential(self, unit: str = 'kcal/mol', as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, expressed in *unit*, and return a dictionary of Numpy Arrays.

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
            return self._get_sigma('SIGMAPOTENTIAL', 'mu', unit=unit, as_df=as_df)
        except KeyError:
            return self._get_sigma('PURESIGMAPOTENTIAL', 'mu', unit=unit, as_df=as_df)

    def get_solubillity(self, subsection: str, as_df: bool = True) -> dict:
        r""""Grab all (temperature dependant) solubillities, returning a dictionary of Numpy Arrays.

        The temperature is stored under the ``"T (K)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self._get_sol('SOLUBILITY', subsection, as_df=as_df)
        except KeyError:
            return self._get_sol('PURESOLUBILITY', subsection, as_df=as_df)

    def get_boiling_point(self, subsection: str, as_df: bool = True) -> dict:
        r""""Grab all (pressure dependant) boiling points, returning a dictionary of Numpy Arrays.

        The pressure is stored under the ``"P (bar)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self._get_bp('BOILINGPOINT', subsection, as_df=as_df)
        except KeyError:
            return self._get_bp('PUREBOILINGPOINT', subsection, as_df=as_df)

    def get_vapor_pressure(self, subsection: str, as_df: bool = True) -> dict:
        r""""Grab all (temperature dependant) solubillities, returning a dictionary of Numpy Arrays.

        The pressure is stored under the ``"T (K)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.
        The latter allows for more convenient plotting of the sigma profile with the
        :meth:`DataFrame.plot` method.

        .. note::
            *as_df* = ``True`` requires Pandas_; data plotting furthermore requires Matplotlib_.

            .. _Pandas: https://pandas.pydata.org/
            .. _Matplotlib: https://matplotlib.org/index.html

        """
        try:
            return self._get_sol('VAPORPRESSURE', subsection, as_df=as_df)
        except KeyError:
            return self._get_sol('PUREVAPORPRESSURE', subsection, as_df=as_df)

    def _get_bp(self, section: str, subsection: str, as_df: bool = False) -> dict:
        """ """
        ret = self.get_array_dict(section, subsection)
        ret['P (bar)'] = self.readarray(section, 'pressure', dtype=float)

        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section, 'P (bar)')

    def _get_sol(self, section: str, subsection: str, as_df: bool = False) -> dict:
        """ """
        ret = self.get_array_dict(section, subsection)
        ret['T (K)'] = self.readarray(section, 'temperature', dtype=float)

        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section, 'T (K)')

    def _get_sigma(self, section: str, subsection: str,
                   unit: str = 'kcal/mol', as_df: bool = False) -> dict:
        """Grab all values of sigma and the sigmapotential/profile;
        combine them into a dictionary or pandas dataframe.

        """
        ret = self.get_array_dict(section, subsection, unit)
        ret['σ (e/A**2)'] = self.readarray(section, 'chdval', dtype=float)

        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section, 'σ (e/A**2)')

    def get_array_dict(self, section: str, subsection: str, unit: str = 'kcal/mol') -> dict:
        """Construct all values for the y-axis"""
        # Use filenames as keys
        _filenames = self.readkf(section, 'filename').split()
        filenames = [_filenames] if not isinstance(_filenames, list) else _filenames
        keys = [os.path.basename(key) for key in filenames]

        # Use sigma profiles/potentials as values
        ratio = Units.conversion_ratio('kcal/mol', unit)
        values = ratio * self.readarray(section, subsection, dtype=float)
        values.shape = len(keys), len(values) // len(keys)

        ret = dict(zip(keys, values))
        if 'PURE' not in section:  # Add a final key for the solvent mixture
            ret['Total'] = self.readarray(section, subsection + 'tot', dtype=float)
        return ret

    @staticmethod
    def _dict_to_df(sigm_dict: dict, section: str, index_name: str) -> 'pandas.DataFrame':
        """Attempt to convert a dictionary, produced by :meth:`._get_sigma` into a DataFrame."""
        if not PANDAS:
            raise ImportError("CRSResults._get_sigma: as_df=True requires the 'pandas' package")

        index = pd.Index(sigm_dict.pop(index_name), name=index_name)
        df = pd.DataFrame(sigm_dict, index=index)
        df.columns.name = section.lower()
        return df


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
