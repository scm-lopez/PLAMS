import os
import inspect
import subprocess

import numpy as np

try:
    import pandas as pd
    PANDAS = True
except ImportError:
    PANDAS = False

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False

from ...tools.units import Units
from .scmjob import (SCMJob, SCMResults)

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

    def get_sigma_profile(self, subsection: str = 'profil', as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, returning a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"σ (e/A**2)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, 'σ (e/A**2)', 'chdval')
        try:
            return self._get_array_dict('SIGMAPROFILE', *args, as_df=as_df)
        except KeyError:
            return self._get_array_dict('PURESIGMAPROFILE', *args, as_df=as_df)

    def get_sigma_potential(self, subsection: str = 'mu', unit: str = 'kcal/mol',
                            as_df: bool = False) -> dict:
        r"""Grab all sigma profiles, expressed in *unit*, and return a dictionary of Numpy Arrays.

        Values of :math:`\sigma` are stored under the ``"σ (e/A**2)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, 'σ (e/A**2)', 'chdval')
        try:
            return self._get_array_dict('SIGMAPOTENTIAL', *args, unit=unit, as_df=as_df)
        except KeyError:
            return self._get_array_dict('PURESIGMAPOTENTIAL', *args, unit=unit, as_df=as_df)

    def get_solubility(self, subsection: str = 'mol_per_L_solvent', as_df: bool = False) -> dict:
        r""""Grab all (temperature dependant) solubilities, returning a dictionary of Numpy Arrays.

        The temperature is stored under the ``"T (K)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, 'T (K)', 'temperature')
        try:
            return self._get_array_dict('SOLUBILITY', *args, as_df=as_df)
        except KeyError:
            return self._get_array_dict('PURESOLUBILITY', *args, as_df=as_df)

    def get_boiling_point(self, subsection: str = 'temperature', as_df: bool = False) -> dict:
        r""""Grab all (pressure dependant) boiling points, returning a dictionary of Numpy Arrays.

        The pressure is stored under the ``"P (bar)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = ('temperature', 'P (bar)', 'pressure')
        try:
            return self._get_array_dict('BOILINGPOINT', *args, as_df=as_df)
        except KeyError:
            return self._get_array_dict('PUREBOILINGPOINT', *args, as_df=as_df)

    def get_vapor_pressure(self, subsection: str = 'vapor pressure', as_df: bool = False) -> dict:
        r""""Grab all (temperature dependant) solubillities, returning a dictionary of Numpy Arrays.

        The pressure is stored under the ``"T (K)"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """
        args = (subsection, 'T (K)', 'temperature')
        try:
            return self._get_array_dict('VAPORPRESSURE', *args, as_df=as_df)
        except KeyError:
            return self._get_array_dict('PUREVAPORPRESSURE', *args, as_df=as_df)

    def get_bi_mixture(self, subsection: str = 'gamma', unit: str = 'kcal/mol',
                       as_df: bool = False) -> dict:
        r""""Grab all (ratio dependant) activity coefficients of a binary mixture, returning a dictionary of Numpy Arrays.

        The component ratio is stored under the ``"molar ratio"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """  # noqa
        args = (subsection, 'molar ratio', 'molar fraction')
        return self._get_array_dict('BINMIXCOEF', *args, unit=unit, as_df=as_df)

    def get_tri_mixture(self, subsection: str = 'gamma', unit: str = 'kcal/mol',
                        as_df: bool = False) -> dict:
        r""""Grab all (ratio dependant) activity coefficients of a ternary mixture, returning a dictionary of Numpy Arrays.

        The component ratio is stored under the ``"molar ratio"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """  # noqa
        args = (subsection, 'molar ratio', 'molar fraction')
        return self._get_array_dict('TERNARYMIX', *args, unit=unit, as_df=as_df)

    def get_composition_line(self, subsection: str = 'gamma', unit: str = 'kcal/mol',
                             as_df: bool = False) -> dict:
        r""""Perform a linear interpolation between the activity coefficients of two components, returning a dictionary of Numpy Arrays.

        The component ratio is stored under the ``"molar ratio"`` key.

        Results can be returned as a Pandas DataFrame by settings *as_df* to ``True``.

        The returned results can be plotted by passing them to the :meth:`CRSResults.plot` method.

        .. note::
            *as_df* = ``True`` requires the Pandas_ package.
            Plotting requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. _Pandas: https://pandas.pydata.org/

        """  # noqa
        args = (subsection, 'molar ratio', 'solvent fraction')
        return self._get_array_dict('COMPOSITIONLINE', *args, unit=unit, as_df=as_df)

    @classmethod
    def plot(cls, array_dict: dict, index_name: str = None, plot_fig: bool = True) -> 'matplotlib.figure.Figure':
        """Plot, show and return a series of COSMO-RS results as a matplotlib Figure instance.

        Accepts the output of, *e.g.*, :meth:`CRSResults.get_sigma_profile`:
        A dictionary of Numpy arrays or a Pandas DataFrame.

        Returns a matplotlib Figure_ instance which can be further modified to the users liking.
        Automatic plotting of the resulting figure can be disabled with the *plot_fig* argument.

        .. note::
            This method requires the `matplotlib <https://matplotlib.org/index.html>`_ package.

        .. note::
            The name of the dictionary/DataFrame key containing the index (*i.e.* the x-axis) can,
            and should, be manually specified in *index_name* if a custom *index_name* is passed
            to :meth:`CRSResults._get_array_dict`.
            This argument can be ignored otherwise.

        .. _Figure: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html#matplotlib.figure.Figure

        """  # noqa
        def get_index_name(array_dict, index_name):
            """Find and return the index and its name."""
            if index_name is not None:
                return array_dict[index_name], index_name

            if not isinstance(array_dict, dict):
                # *array_dict* is a dataframe instead of a dictionary
                return array_dict.index, array_dict.index.name

            idx_tups = ('σ (e/A**2)', 'T (K)', 'P (bar)', 'molar ratio')
            for k, v in array_dict.items():
                if k in idx_tups:
                    return v, k

        # Check if matplotlib is installed
        if not MATPLOTLIB:
            method = cls.__name__ + '.plot'
            raise ImportError("{}: this method requires the 'matplotlib' package".format(method))

        # Retrieve the index and its name
        index, index_name = get_index_name(array_dict, index_name)

        # Assign various series to the plot
        fig, ax = plt.subplots()
        for k, v in array_dict.items():
            if k == index_name:
                continue
            ax.plot(index, v, label=k)

        # Add the legend and x-label
        ax.legend()
        ax.set_xlabel(k)
        try:
            ax.set_title(array_dict.columns.name)
        except AttributeError:
            pass  # *array_dict* is not a dataframe

        # Show and return
        if plot_fig:
            plt.show()
        return fig

    def _get_array_dict(self, section: str, subsection: str, index_name: str, index_subsection: str,
                        unit: str = 'kcal/mol', as_df: bool = False) -> dict:
        """Create dictionary or DataFrame containing all values in *section*/*subsection*.

        Takes the following arguments:
            * The *section*/*subsection* of the desired quantity.
            * The desired name of the index (*index_name*).
            * The name of subsection containing the index (*index_subsection*).
            * The *unit* of the output quanty (ignore this keyword if not applicable).
            * If the result should be returned as Pandas DataFrame (*as_df*).

        """
        ret = self._construct_array_dict(section, subsection, unit)

        # Create the index
        index = self.readarray(section, index_subsection, dtype=float)
        if section in ('BINMIXCOEF', 'COMPOSITIONLINE', 'TERNARYMIX'):
            ncomponent = 3 if section == 'TERNARYMIX' else 2
            index.shape = ncomponent, len(index) // ncomponent
            iterator = np.nditer(index.astype(str), flags=['external_loop'], order='F')
            ret[index_name] = np.array([' / '.join(i for i in item) for item in iterator])
        else:
            ret[index_name] = index

        # Return a dictionary of arrays or a DataFrame
        if not as_df:
            return ret
        else:
            return self._dict_to_df(ret, section, index_name)

    def _construct_array_dict(self, section: str, subsection: str, unit: str = 'kcal/mol') -> dict:
        """Construct dictionary containing all values in *section*/*subsection*."""
        # Use filenames as keys
        _filenames = self.readkf(section, 'filename').split()
        filenames = [_filenames] if not isinstance(_filenames, list) else _filenames

        # Grab the keys and the number of items per key
        keys = [os.path.basename(key) for key in filenames] + ['Total']
        nitems = self.readkf(section, 'nitems')

        # Use sigma profiles/potentials as values
        ratio = Units.conversion_ratio('kcal/mol', unit)
        values = ratio * self.readarray(section, subsection, dtype=float)
        values.shape = len(values) // nitems, nitems

        ret = dict(zip(keys, values))
        try:
            ret['Total'] = self.readarray(section, subsection + 'tot', dtype=float)
        except KeyError:
            pass
        return ret

    @classmethod
    def _dict_to_df(cls, array_dict: dict, section: str, index_name: str) -> 'pandas.DataFrame':
        """Attempt to convert a dictionary into a DataFrame."""
        if not PANDAS:
            method = cls.__name__ + '.' + inspect.stack()[2][3]
            raise ImportError("{}: as_df=True requires the 'pandas' package".format(method))

        index = pd.Index(array_dict.pop(index_name), name=index_name)
        df = pd.DataFrame(array_dict, index=index)
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

    @classmethod
    def cos_to_coskf(cls, filename: str) -> str:
        """Convert a .cos file into a .coskf file with the :code:`$ADFBIN/cosmo2kf` command.

        Returns the filename of the new .coskf file.

        """
        filename_out = filename + 'kf'
        try:
            adfbin = os.environ['ADFBIN']
        except KeyError:
            method = cls.__name__ + '.cos_to_coskf'
            err = '{}: the "ADFBIN" environment variable was not found'.format(method)
            raise EnvironmentError(err)

        args = [os.path.join(adfbin, 'cosmo2kf'), filename, filename_out]
        subprocess.run(args)
        return filename_out
