"""A module with various context managers."""

import threading
from abc import ABC, abstractmethod
from functools import wraps
from typing import Any, ContextManager, Collection, Callable, Dict

from .errors import ReentranceError

__all__ = []


class FuncReplacerABC(ContextManager[None], ABC):
    """An abstract context manager for temporary replacing one or more methods of a class.

    Subclasses will have to defined two attributes:

    * :attr:`FuncReplacerABC.replace_func`: A :class:`~collections.abc.Collection` with the
      names of all to-be replaced methods.
    * :meth:`FuncReplacerABC.decorate`: A static method for altering/replacing the methods.

    Instances of this class are thread-safe (note that this does *not* include
    :meth:`FuncReplacerABC.__init__` itself), resuable, but non-reentrant:

    .. code:: python

        >>> class Context(FuncReplacerABC):
        ...     ...

        >>> manager = Context()

        # Opening the context manager multiple times is completelly fine
        >>> with manager:
        ...     ...
        >>> with manager:
        ...     ...
        >>> with manager:
        ...     ...

        # This is not fine
        >>> with manager:
            >>> with manager:
            ...     ...
        ReentranceError("'Context' instances cannot not be entered in a reentrant manner")

    """

    #: A private instance variable for ensuring the thread safety of
    #: :meth:`__enter__` and :meth:`__exit__` calls.
    _lock: threading.Lock

    #: A private instance variable for keeping track of (potentially reentrant)
    #: :meth:`__enter__` calls.
    _open: bool

    #: An instance variable containing object whose methods are to-be replaced upon
    #: calling :meth:`__enter__`.
    obj: type

    #: An instance variable containing a dictionary mapping the names of the to-be replaced
    #: functions to the old functions.
    func_old: Dict[str, Callable]

    #: An instance variable containing a dictionary mapping the names of the to-be replaced
    #: functions to the new functions.
    func_new: Dict[str, Callable]

    @property
    @abstractmethod
    def replace_func(self) -> Collection[str]:
        """A class variable defining all to-be replaced attributes for when the context manager has been entered.

        This attribute should be defined in a subclass as following:

        .. code:: python

            >>> class SubClass(FuncReplacerABC):
            ...     replace_func = ('name1', 'name2', 'name3', ...)  # Any Collection is fine

        """
        raise NotImplementedError('Trying to get an abstract attribute')

    @staticmethod
    @abstractmethod
    def decorate(func: Callable) -> Callable:
        """Decorate all functions defined in :attr:`replace_func`."""
        raise NotImplementedError('Trying to call an abstract attribute')

    def __init__(self, obj: type) -> None:
        """Initialize the context manager.

        Note that, contrary to :meth:`__enter__` and :meth:`__exit__`,
        this method is *not* thread-safe.

        """
        self._lock = threading.Lock()
        self._open = False

        # Ensure that obj is a class, not a class instance
        type_obj = obj if isinstance(obj, type) else type(obj)
        self.obj = type_obj

        # Define the old and new method
        self.func_old = {name: getattr(type_obj, name) for name in self.replace_func}
        self.func_new = {name: self.decorate(func) for name, func in self.func_old.items()}

    def __eq__(self, obj: Any) -> bool:
        """Implement :code:`self == obj`."""
        if type(self) is not type(obj):
            return False
        return self.obj is getattr(obj, 'obj', None)

    def __repr__(self) -> str:
        """Implement :code:`repr(self)` and :code:`str(self)`."""
        return f'{self.__class__.__name__}(obj={self.obj!r})'

    def __enter__(self) -> None:
        """Enter the context manager: modify all methods in :attr:`func_new` at the class level."""
        # Precaution against calling __enter__() in a recursive manner
        with self._lock:
            if self._open:
                raise ReentranceError(f"{self.__class__.__name__!r} instances cannot "
                                      "not be entered in a reentrant manner")

            for name, func in self.func_new.items():
                setattr(self.obj, name, func)
            self._open = False

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager: restore all methods in :attr:`func_old` at the class level."""
        if exc_type is ReentranceError:
            return

        with self._lock:
            for name, func in self.func_old.items():
                setattr(self.obj, name, func)
            self._open = False


class SupressMissing(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary disabling the :meth:`__missing__` magic method.

    See :meth:`Settings.supress_missing` for more details.

    """  # noqa

    replace_func = ('__missing__',)

    @staticmethod
    def decorate(func):
        """Decorate *func* such that calling it raises a :exc:`KeyError`."""
        @wraps(func)
        def wrapper(self, name):
            raise KeyError(name)
        return wrapper


class Lower(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary converting all keys passed to :meth:`__delitem__`, :meth:`__setitem__` and :meth:`__getitem__` to lower case."""

    replace_func = ('__delitem__', '__setitem__', '__getitem__')

    @staticmethod
    def decorate(func):
        """Decorate *func* such that the passed key is converted to lower case before being called."""
        @wraps(func)
        def wrapper(self, key, *args, **kwargs):
            try:
                key_low = key.lower()
            except (AttributeError, TypeError):
                key_low = key
            return func(self, key_low, *args, **kwargs)
        return wrapper


class Upper(FuncReplacerABC):
    """A reusable, but non-reentrant, context manager for temporary converting all keys passed to :meth:`__delitem__`, :meth:`__setitem__` and :meth:`__getitem__` to upper case."""

    replace_func = ('__delitem__', '__setitem__', '__getitem__')

    @staticmethod
    def decorate(func):
        """Decorate *func* such that the passed key is converted to upper case before being called."""
        @wraps(func)
        def wrapper(self, key, *args, **kwargs):
            try:
                key_upper = key.upper()
            except (AttributeError, TypeError):
                key_upper = key
            return func(self, key_upper, *args, **kwargs)
        return wrapper
