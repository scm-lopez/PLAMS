import sys
from typing import TYPE_CHECKING, Union, Sequence, Tuple

import numpy as np

# Typing related imports
# Note that all code in the `in TYPE_CHECKING` block will remain unexecuted during runtime
if TYPE_CHECKING:
    from numpy.typing import ArrayLike, DtypeLike

    if sys.version_info[1] >= 8:
        from typing import TypedDict, Literal, SupportsIndex
    else:
        from typing_extensions import TypedDict, Literal, Protocol

        class SupportsIndex(Protocol):
            def __index__(self) -> int:
                pass

    _IndexLike = Union[None, SupportsIndex, Sequence[int], slice, ArrayLike]

    _AttrNames = Literal['atoms', 'bonds', 'mol', 'scale']

    _GetTuple = Tuple[
        bool,  # IsHeteroAtom
        int,  # SerialNumber
        str,  # Name
        str,  # AltLoc
        str,  # ResidueName
        str,  # ChainId
        int,  # ResidueNumber
        str,  # InsertionCode
        float,  # x
        float,  # y
        float,  # z
        float,  # Occupancy
        float,  # TempFactor
        str,  # symbol
        int,  # charge
    ]

    class _PDBInfo(TypedDict):
        IsHeteroAtom: bool
        SerialNumber: int
        Name: str
        AltLoc: str
        ResidueName: str
        ChainId: str
        ResidueNumber: int
        InsertionCode: str
        Occupancy: float
        TempFactor: float

    class _Properties(TypedDict):
        charge: int
        pdb_info: _PDBInfo

    _SetTuple = Tuple[
        _Properties,  # properties
        Tuple[float, float, float],  # coords
        str  # symbol
    ]

else:
    ArrayLike = 'numpy.typing.ArrayLike'
    DtypeLike = 'numpy.typing.DtypeLike'
    SupportsIndex = 'typing.SupportsIndex'

    _IndexLike = f'{__name__}._IndexLike'
    _AttrNames = f'{__name__}._AttrNames'
    _GetTuple = f'{__name__}._GetTuple'
    _SetTuple = f'{__name__}._SetTuple'
