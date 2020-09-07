from typing import (
    TypeVar,
    overload,
    Any,
    ClassVar,
    Container,
    Dict,
    Iterable,
    Iterator,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    Callable,
)

import h5py
import numpy as np

from ._annotations import ArrayLike, DtypeLike, _AttrNames, _IndexLike, Literal
from ...mol.molecule import Molecule

_ST = TypeVar("_ST", bound=HDF5Mol)

class InfoTuple(NamedTuple):
    name: str
    dtype: DtypeLike
    vlen: bool = ...
    ndim: int = ...
    order: Optional[Literal['C', 'F']] = ...

class _FalseType:
    def __eq__(self, b: object) -> Literal[False]: ...

ATOMS_DTYPE: np.dtype
BONDS_DTYPE: np.dtype
MOL_DTYPE: np.dtype
SCALE_DTYPE: np.dtype

class HDF5Mol:
    INFO: ClassVar[Mapping[_AttrNames, InfoTuple]] = ...

    @property
    def atoms(self) -> np.ndarray: ...
    @atoms.setter
    def atoms(self, value: Optional[ArrayLike]) -> None: ...
    @property
    def bonds(self) -> np.ndarray: ...
    @bonds.setter
    def bonds(self, value: Optional[ArrayLike]) -> None: ...
    @property
    def mol(self) -> np.ndarray: ...
    @mol.setter
    def mol(self, value: Optional[ArrayLike]) -> None: ...
    @property
    def scale(self) -> np.ndarray: ...
    @scale.setter
    def scale(self, value: Optional[ArrayLike]) -> None: ...

    def __init_subclass__(cls) -> None: ...
    @overload
    def __init__(
        self,
        atoms: np.ndarray,
        bonds: np.ndarray,
        mol: np.ndarray,
        scale: np.ndarray,
        validate: Literal[False],
        copy: bool = ...
    ) -> None: ...
    @overload
    def __init__(
        self,
        atoms: Optional[ArrayLike] = ...,
        bonds: Optional[ArrayLike] = ...,
        mol: Optional[ArrayLike] = ...,
        scale: Optional[ArrayLike] = ...,
        validate: Literal[True] = ...,
        copy: bool = ...,
    ) -> None: ...
    @overload
    def _init(
        self,
        array_dict: Mapping[_AttrNames, np.ndarray],
        validate: Literal[False],
        copy: bool = ...
    ) -> None: ...
    @overload
    def _init(
        self,
        array_dict: Mapping[_AttrNames, Optional[ArrayLike]],
        validate: Literal[True] = ...,
        copy: bool = ...,
    ) -> None: ...
    def __reduce__(
        self,
    ) -> Tuple[Type[_ST], Tuple[Literal[False]], Dict[_AttrNames, np.ndarray]]: ...
    def __setstate__(self, state: Dict[_AttrNames, np.ndarray]) -> None: ...
    def __copy__(self) -> _ST: ...
    def __deepcopy__(self, memo: Any = ...) -> _ST: ...
    def __len__(self) -> int: ...
    def __eq__(self, value: object) -> bool: ...
    def __getitem__(self, index: _IndexLike) -> _ST: ...

    def copy(self) -> _ST: ...
    @overload
    @staticmethod
    def _eq(vlen: Literal[True], ar1: np.ndarray, ar2: np.ndarray) -> bool: ...
    @overload
    @staticmethod
    def _eq(vlen: Literal[False], ar1: np.ndarray, ar2: np.ndarray) -> np.bool_: ...
    @classmethod
    def _parse_array(
        cls, name: _AttrNames, value: Optional[ArrayLike], size: int, copy: bool = ...
    ) -> np.ndarray: ...
    @classmethod
    def _parse_scale(
        cls, name: _AttrNames, value: Optional[ArrayLike], size: int, copy: bool = ...
    ) -> np.ndarray: ...

    @classmethod
    def keys(cls) -> Iterator[_AttrNames]: ...
    def values(self) -> Iterator[np.ndarray]: ...
    def items(self) -> Iterator[Tuple[_AttrNames, np.ndarray]]: ...
    def concatenate(self, *args: _ST) -> _ST: ...
    @classmethod
    def from_molecules(
        cls: Type[_ST],
        mol_list: Sequence[Molecule],
        scale: Optional[ArrayLike] = ...,
        keys: Optional[Container[_AttrNames]] = ...,
    ) -> _ST: ...
    @classmethod
    def _from_molecules(
        cls: Type[_ST],
        mol_list: Sequence[Molecule],
        array_dict: Mapping[_AttrNames, Optional[ArrayLike]],
        keys: Optional[Container[_AttrNames]] = ...,
    ) -> _ST: ...
    @classmethod
    def _get_info(
        cls,
        mol_list: Sequence[Molecule],
        name: _AttrNames,
        func: Callable[[Molecule, DtypeLike], np.ndarray],
    ) -> np.ndarray: ...
    @staticmethod
    def get_atom_info(mol: Molecule, dtype: DtypeLike) -> np.ndarray: ...
    @staticmethod
    def get_bond_info(mol: Molecule, dtype: DtypeLike) -> np.ndarray: ...
    @staticmethod
    def get_mol_info(mol: Molecule, dtype: DtypeLike) -> np.ndarray: ...
    def to_molecules(
        self,
        mol: Union[None, Molecule, Iterable[Optional[Molecule]]] = ...,
        keys: Optional[Container[_AttrNames]] = ...,
    ) -> List[Molecule]: ...
    @classmethod
    def _to_molecules(
        cls,
        atoms: np.ndarray,
        bonds: np.ndarray,
        mol: np.ndarray,
        keys: Container[_AttrNames],
        molecule: Optional[Molecule] = ...,
    ) -> Molecule: ...
    @staticmethod
    def set_atom_info(mol: Molecule, data: np.ndarray) -> None: ...
    @staticmethod
    def set_bond_info(mol: Molecule, data: np.ndarray) -> None: ...
    @staticmethod
    def set_mol_info(mol: Molecule, data: np.void) -> None: ...
    def create_group(
        self,
        file: h5py.Group,
        name: str,
        atoms: Optional[h5py.Dataset] = ...,
        bonds: Optional[h5py.Dataset] = ...,
        mol: Optional[h5py.Dataset] = ...,
        scale: Optional[h5py.Dataset] = ...,
        create_group: bool = ...,
        **kwargs: Any
    ) -> h5py.Group: ...
    def _create_group(
        self,
        file: h5py.Group,
        name: str,
        dset_dict: Mapping[_AttrNames, Optional[h5py.Dataset]],
        create_group: bool = ...,
        **kwargs: Any
    ) -> h5py.Group: ...
    def to_hdf5(
        self,
        group: h5py.Group,
        index: _IndexLike = ...,
        keys: Optional[Container[_AttrNames]] = ...,
    ) -> None: ...
    def _update_hdf5_shape(self, group: h5py.Group, min_len: int) -> None: ...
    @classmethod
    def from_hdf5(
        cls: Type[_ST],
        group: h5py.Group,
        index: _IndexLike = ...,
        keys: Optional[Container[_AttrNames]] = ...,
    ) -> _ST: ...
    def intersection(self, value: Union[_ST, ArrayLike]) -> _ST: ...
    def difference(self, value: Union[_ST, ArrayLike]) -> _ST: ...
    def symmetric_difference(self, value: _ST) -> _ST: ...
    def union(self, value: _ST) -> _ST: ...
    def _get_index(
        self: _ST, value: Union[_ST, ArrayLike]
    ) -> Tuple[np.ndarray, np.ndarray]: ...
