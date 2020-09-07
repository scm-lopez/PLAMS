import types
import inspect
import functools
from typing import Union, overload, Optional, Callable, Iterable

import numpy as np

from ._annotations import ArrayLike, DtypeLike, SupportsIndex

try:
    import h5py
except ImportError as ex:
    H5PY_EX: Optional[ImportError] = ex
else:
    H5PY_EX = None


def _alen(a: ArrayLike) -> int:
    """Backport of the deprecated :func:`numpy.alen` function."""
    try:
        return len(a)
    except TypeError:
        arr = np.atleast_1d(a)
        return len(arr)


@overload
def _parse_index(index: Union[None, slice, SupportsIndex], seq_len: Optional[int]) -> slice:
    ...
@overload
def _parse_index(index: ArrayLike, seq_len: Optional[int]) -> np.ndarray:
    ...
def _parse_index(index, seq_len):
    """Convert the passed object into a :class:`slice` or :class:`numpy.ndarray`."""
    if index is None:
        return slice(None)
    elif isinstance(index, slice):
        return index

    try:
        return _int_to_slice(index, seq_len)
    except (AttributeError, TypeError):
        ret = np.atleast_1d(index)
        if ret.ndim != 1:
            raise ValueError("'index' expected an array with 1 dimension; "
                             f"observed dimensionality: {ret.ndim}") from None
        return ret


def _int_to_slice(int_like: SupportsIndex, seq_len: int) -> slice:
    """Take an integer-like object and convert it into a :class:`slice`."""
    integer = int_like.__index__()
    if integer > 0:
        if integer != seq_len:
            return slice(None, integer + 1)
        else:
            return slice(integer - 1, None)

    else:
        if integer == -1:
            return slice(integer, None)
        else:
            return slice(integer, integer + 1)


def _to_array(value: Optional[ArrayLike], size: int, dtype: DtypeLike,
              vlen: bool = False,
              copy: bool = False,
              ndim: int = 1,
              order: Optional[str] = None) -> np.ndarray:
    """Convert **value** into an array of the specified **size** and **dtype**."""
    if H5PY_EX is not None:
        raise H5PY_EX

    # Create and return an empty array
    if value is None:
        if vlen:
            shape = (0,) + (1,) * (ndim - 1)
            vlen_dtype = h5py.vlen_dtype(dtype)
            array = np.empty(size, dtype=vlen_dtype)
            array[:] = size * [np.zeros(shape, dtype=dtype, order=order)]
            return array
        else:
            shape = (size,) + (1,) * (ndim - 1)
            return np.zeros(shape, dtype=dtype, order=order)

    # Create a ragged array if `vlen is True`
    if not vlen:
        ret = np.array(value, dtype=dtype, copy=copy, order=order)

        # Broadcast the imput if `value` is a scalar
        if not ret.ndim:
            return np.full(shape, ret, dtype=dtype, order=order)
        ret.shape += (1,) * (ndim - ret.ndim)

    else:
        vlen_dtype = h5py.vlen_dtype(dtype)
        try:
            return np.asarray(value, vlen_dtype)
        except ValueError:
            raise TypeError("The setting of variable-length arrays is not supported")

    if len(ret) != size:
        raise ValueError(f'The passed array should be of length {size}')
    return ret
