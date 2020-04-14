"""A module with various context managers."""

import threading
from abc import abstractmethod, ABCMeta
from types import MappingProxyType
from weakref import WeakValueDictionary
from functools import wraps
from typing import (
    Any, Collection, Callable, Dict, Tuple, NoReturn, ContextManager, cast,
    Type, Optional, ClassVar, MutableMapping, Mapping, TypeVar, Generic, Union
)

from .errors import ReentranceError

__all__ = ['FuncReplacerABC']

FT = TypeVar('FT', bound='FuncReplacerABC')


class _FuncReplacerMeta(ABCMeta):
    """The metaclass of :class:`FuncReplacerABC`.

    Used for setting the :attr:`FuncReplacerABC._type_cache` class variable,
    a :class:`~weakref.WeakValueDictionary` for keeping track of all class instance singletons.

    """

    def __new__(mcls, name, bases, namespace, **kwargs):
        cls = super().__new__(mcls, name, bases, namespace, **kwargs)
        cls._type_cache = WeakValueDictionary()
        return cls


class FuncReplacerABC(ContextManager[None], metaclass=_FuncReplacerMeta):
    """An abstract context manager for temporary replacing one or more methods of a class.

    Subclasses will have to defined two attributes:

    * :attr:`FuncReplacerABC.replace_func`: A :class:`~collections.abc.Collection` with the
      names of all to-be replaced methods.
    * :meth:`FuncReplacerABC.decorate`: A decorator for altering/replacing the methods.

    Instances of this class furthermore have the following properties:

    * Entering and exiting is thread-safe (note that this does *not* include
      :meth:`FuncReplacerABC.__init__` itself).
    * They are immutable(-ish).
    * They are singletons with respect to the passed object type.
    * They are resuable but non-reentrant.

    For example:

    .. code:: python

        # Define and instantiate a subclass
        >>> class Context(FuncReplacerABC):
        ...     replace_func = ('__int__', ...)
        ...
        ...     @staticmethod
        ...     def decorate(func): return ...

        # Opening the context manager multiple times is completely fine
        >>> with manager:
        ...     ...
        >>> with manager:
        ...     ...
        >>> with manager:
        ...     ...

        # This is not fine
        >>> with manager:
        ...     with manager:
        ...         ...
        ReentranceError: "'Context' instances cannot not be entered in a reentrant manner"

        # Class instances are singletons with respect to the passed object type
        >>> manager1 = Context(int)
        >>> manager2 = Context(1)
        >>> manager3 = Context(999)
        >>> print(manager1 is manager2 is manager3)
        True

        # Fiddling with attributes is not fine
        >>> del manager.obj
        AttributeError: attribute 'obj' of 'Context' objects is not writable

    """
    #: A class variable for keeping track of all :class:`FuncReplacerABC` instances.
    #: Used for ensuring all instances are singletons with respect to the passed object type.
    _type_cache: ClassVar[MutableMapping[type, 'FuncReplacerABC']]

    #: An instance variable for ensuring the thread safety of
    #: :meth:`__enter__` and :meth:`__exit__` calls.
    _lock: threading.Lock

    #: An instance variable for keeping track of (potentially reentrant)
    #: :meth:`__enter__` calls.
    _open: bool

    #: An instance variable containing object whose methods are to-be replaced upon
    #: calling :meth:`__enter__`.
    obj: type

    #: An instance variable containing a read-only :class:`~collections.abc.Mapping` whichs
    #: maps the names of the to-be replaced functions to the old functions.
    func_old: Mapping[str, Callable]

    #: An instance variable containing a read-only :class:`~collections.abc.Mapping` whichs
    #: maps the names of the to-be replaced functions to the new functions.
    func_new: Mapping[str, Callable]

    # Abstract methods and attributes

    @property
    @abstractmethod
    def replace_func(self) -> Collection[str]:
        """A class variable defining all to-be replaced attributes for when the context manager has been entered.

        This attribute should be defined in a subclass as following:

        .. code:: python

            >>> class SubClass(FuncReplacerABC):
            ...     replace_func = ('name1', 'name2', 'name3', ...)  # Any Collection is fine

        """  # noqa: E501
        raise NotImplementedError('Trying to get an abstract attribute')

    @staticmethod
    @abstractmethod
    def decorate(func: Callable) -> Callable:
        """Decorate all functions defined in :attr:`replace_func`."""
        raise NotImplementedError('Trying to call an abstract method')

    # Various magic methods

    def __new__(cls: Type[FT], obj: Union[type, Any]) -> FT:
        """Construct a new :class:`FuncReplacerABC` instance."""
        # If possible, return a cached FuncReplacerABC instance
        type_obj = cast(type, obj if isinstance(obj, type) else type(obj))
        try:
            return cast(FT, cls._type_cache[type_obj])
        except KeyError:
            return super().__new__(cls)

    def __init__(self, obj: Union[type, Any]) -> None:
        """Initialize the context manager.

        Note that, contrary to :meth:`__enter__` and :meth:`__exit__`,
        this method is *not* thread-safe.

        """
        # Check if self is a cached (already initialized) instance
        if vars(self):
            return
        cls = type(self)

        setattr = super().__setattr__
        setattr('_lock', threading.Lock())
        setattr('_open', False)

        # Ensure that obj is a class, not a class instance
        type_obj = cast(type, obj if isinstance(obj, type) else type(obj))
        setattr('obj', type_obj)

        # Define the old and new method
        setattr('func_old', MappingProxyType({name: getattr(type_obj, name) for
                                              name in cls.replace_func}))  # type: ignore
        setattr('func_new', MappingProxyType({name: self.decorate(func) for
                                              name, func in self.func_old.items()}))

        # Update the type cache
        cls._type_cache[type_obj] = self

    def __reduce__(self: FT) -> Tuple[Type[FT], Tuple[type]]:
        """Helper function for :mod:`pickle`."""
        return (type(self), (self.obj,))

    def __copy__(self: FT) -> FT:
        """Implement :code:`copy.copy`."""
        return self

    def __deepcopy__(self: FT, memo: Optional[Dict[int, Any]] = None) -> FT:
        """Implement :code:`copy.deepcopy`."""
        return self

    def __repr__(self) -> str:
        """Implement :code:`repr(self)` and :code:`str(self)`."""
        return f'{self.__class__.__name__}(obj={self.obj!r})'

    def __setattr__(self, name: str, value: Any) -> NoReturn:
        """Implement :code:`setattr(self, name, value)`."""
        raise self._attributeError(name)  # Attributes are read-only

    def __delattr__(self, name: str) -> NoReturn:
        """Implement :code:`delattr(self, name)`."""
        raise self._attributeError(name)  # Attributes are read-only

    def _attributeError(self, name: str) -> AttributeError:
        """Return an :exc:`AttributeError`; attributes of this instance are read-only."""
        if hasattr(self, name):
            cls_name = self.__class__.__name__
            return AttributeError(f"attribute {name!r} of {cls_name!r} objects is not writable")
        else:
            return AttributeError(f"{self.__class__.__name__!r} object has no attribute {name!r}")

    # Context manager-related magic methods

    def __enter__(self) -> None:
        """Enter the context manager: modify all methods in :attr:`func_new` at the class level."""
        if self._open:
            raise ReentranceError(f"{self.__class__.__name__!r} instances cannot "
                                  "not be entered in a reentrant manner")

        # Precaution against calling __enter__() in a recursive manner
        with self._lock:
            for name, func in self.func_new.items():
                setattr(self.obj, name, func)
        super().__setattr__('_open', True)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager: restore all methods in :attr:`func_old` at the class level."""
        if exc_type is ReentranceError:
            return

        with self._lock:
            for name, func in self.func_old.items():
                setattr(self.obj, name, func)
        super().__setattr__('_open', False)


class SupressMissing(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary disabling the :meth:`__missing__` magic method.

    See :meth:`Settings.supress_missing` for more details.

    """  # noqa: E501

    replace_func = ('__missing__',)

    @staticmethod
    def decorate(func):
        """Decorate *func* such that calling it raises a :exc:`KeyError`."""
        @wraps(func)
        def wrapper(self, name):
            raise KeyError(name)
        return wrapper


# Work in progress context managers


class _Lower(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary converting all keys to lower case."""  # noqa: E501

    replace_func = ('__delitem__', '__setitem__', '__getitem__', '__contains__',
                    'get', 'pop', 'popitem', 'setdefault')

    @staticmethod
    def decorate(func):
        """Decorate *func* such that the passed key is converted to lower case before being called."""  # noqa: E501
        @wraps(func)
        def wrapper(self, key, *args, **kwargs):
            try:
                key_low = key.lower()
            except (AttributeError, TypeError):
                key_low = key
            return func(self, key_low, *args, **kwargs)
        return wrapper


class _Upper(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary converting all keys to upper case."""  # noqa: E501

    replace_func = ('__delitem__', '__setitem__', '__getitem__', '__contains__',
                    'get', 'pop', 'popitem', 'setdefault')

    @staticmethod
    def decorate(func):
        """Decorate *func* such that the passed key is converted to upper case before being called."""  # noqa: E501
        @wraps(func)
        def wrapper(self, key, *args, **kwargs):
            try:
                key_upper = key.upper()
            except (AttributeError, TypeError):
                key_upper = key
            return func(self, key_upper, *args, **kwargs)
        return wrapper
