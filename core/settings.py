import textwrap
import contextlib

from .context import SupressMissing, Lower, Upper

__all__ = ['Settings', 'ig']


class _SettingsMeta(type):
    def __new__(mcls, name, bases, namespace, **kwargs):
        cls = super().__new__(mcls, name, bases, namespace, **kwargs)
        cls._supress_missing = SupressMissing(cls)
        cls._lower = Lower(cls)
        cls._upper = Upper(cls)
        return cls


class Settings(dict, metaclass=_SettingsMeta):
    """Automatic multi-level dictionary. Subclass of built-in class :class:`dict`.

    The shortcut dot notation (``s.basis`` instead of ``s['basis']``) can be used for keys that:

    *   are strings
    *   don't contain whitespaces
    *   begin with a letter or an underscore
    *   don't both begin and end with two or more underscores.

    Iteration follows lexicographical order (via :func:`sorted` function)

    Methods for displaying content (:meth:`~object.__str__` and :meth:`~object.__repr__`) are overridden to recursively show nested instances in easy-readable format.

    Regular dictionaries (also multi-level ones) used as values (or passed to the constructor) are automatically transformed to |Settings| instances::

        >>> s = Settings({'a': {1: 'a1', 2: 'a2'}, 'b': {1: 'b1', 2: 'b2'}})
        >>> s.a[3] = {'x': {12: 'q', 34: 'w'}, 'y': 7}
        >>> print(s)
        a:
          1:    a1
          2:    a2
          3:
            x:
              12:   q
              34:   w
            y:  7
        b:
          1:    b1
          2:    b2

    """

    # To-be replaced with SupressMissing, Lower and Upper instances by the metaclass
    _supress_missing: SupressMissing
    _lower: Lower
    _upper: Upper

    def __init__(self, *args, **kwargs):
        cls = type(self)
        dict.__init__(self, *args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                self[k] = cls(v)
            if isinstance(v, list):
                self[k] = [cls(i) if isinstance(i, dict) else i for i in v]


    def copy(self):
        """Return a new instance that is a copy of this one. Nested |Settings| instances are copied recursively, not linked.

        In practice this method works as a shallow copy: all "proper values" (leaf nodes) in the returned copy point to the same objects as the original instance (unless they are immutable, like ``int`` or ``tuple``). However, nested |Settings| instances (internal nodes) are copied in a deep-copy fashion. In other words, copying a |Settings| instance creates a brand new "tree skeleton" and populates its leaf nodes with values taken directly from the original instance.

        This behavior is illustrated by the following example::

            >>> s = Settings()
            >>> s.a = 'string'
            >>> s.b = ['l','i','s','t']
            >>> s.x.y = 12
            >>> s.x.z = {'s','e','t'}
            >>> c = s.copy()
            >>> s.a += 'word'
            >>> s.b += [3]
            >>> s.x.u = 'new'
            >>> s.x.y += 10
            >>> s.x.z.add(1)
            >>> print(c)
            a:  string
            b:  ['l', 'i', 's', 't', 3]
            x:
              y:    12
              z:    set([1, 's', 'e', 't'])
            >>> print(s)
            a:  stringword
            b:  ['l', 'i', 's', 't', 3]
            x:
              u:    new
              y:    22
              z:    set([1, 's', 'e', 't'])

        This method is also used when :func:`python3:copy.copy` is called.
        """
        cls = type(self)
        ret = cls()
        for name in self:
            if isinstance(self[name], cls):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret



    def soft_update(self, other):
        """Update this instance with data from *other*, but do not overwrite existing keys. Nested |Settings| instances are soft-updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.soft_update(o)
            >>> print(s)
            a:  AA        #original value s.a not overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   XY1   #original value s.x.y1 not overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.

        Shortcut ``A += B`` can be used instead of ``A.soft_update(B)``.
        """
        cls = type(self)
        for name in other:
            if isinstance(other[name], cls):
                if name not in self:
                    self[name] = other[name].copy()
                elif isinstance(self[name], cls):
                    self[name].soft_update(other[name])
            elif name not in self:
                self[name] = other[name]
        return self



    def update(self, other):
        """Update this instance with data from *other*, overwriting existing keys. Nested |Settings| instances are updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.update(o)
            >>> print(s)
            a:  O_AA        #original value s.a overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   O_XY1   #original value s.x.y1 overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.
        """
        cls = type(self)
        for name in other:
            if isinstance(other[name], cls):
                if name not in self or not isinstance(self[name], cls):
                    self[name] = other[name].copy()
                else:
                    self[name].update(other[name])
            else:
                self[name] = other[name]



    def merge(self, other):
        """Return new instance of |Settings| that is a copy of this instance soft-updated with *other*.

        Shortcut ``A + B`` can be used instead of ``A.merge(B)``.
        """
        ret = self.copy()
        ret.soft_update(other)
        return ret



    def find_case(self, key):
        """Check if this instance contains a key consisting of the same letters as *key*, but possibly with different case. If found, return such a key. If not, return *key*.

        When |Settings| are used in case-insensitive contexts, this helps preventing multiple occurences of the same key with different case::

            >>> s = Settings()
            >>> s.system.key1 = 'value1'
            >>> s.System.key2 = 'value2'
            >>> print(s)
            System:
                key2:    value2
            system:
                key1:    value1

            >>> t = Settings()
            >>> t.system.key1 = 'value1'
            >>> t[t.find_case('System')].key2 = 'value2'
            >>> print(t)
            system:
                key1:    value1
                key2:    value2

        """
        lowkey = key.lower()
        for k in self:
            if k.lower() == lowkey:
                return k
        return key



    def as_dict(self):
        """Return a copy of this instance with all |Settings| replaced by regular Python dictionaries.
        """
        cls = type(self)
        d = {}
        for k, v in self.items():
            if isinstance(v, cls):
                d[k] = v.as_dict()
            elif isinstance(v, list):
                d[k] = [i.as_dict() if isinstance(i, cls) else i for i in v]
            else:
                d[k] = v

        return d


    @classmethod
    def supress_missing(cls):
        """A reusable, but non-reentrant, context manager for temporary disabling the :meth:`.Settings.__missing__` magic method: all calls now raising a :exc:`KeyError`.

        As a results, attempting to access keys absent from an arbitrary |Settings| instance will raise a :exc:`KeyError`, thus reverting to the default dictionary behaviour.

        .. note::
            The :meth:`.Settings.__missing__` method is (temporary) supressed at the class level to ensure consistent invokation by the Python interpreter.
            See also `special method lookup`_.

        Example:

        .. code:: python

            >>> s = Settings()

            >>> with s.supress_missing():
            ...     s.a.b.c = True
            KeyError: 'a'

            >>> s.a.b.c = True
            >>> print(s.a.b.c)
            True

        .. _`special method lookup`: https://docs.python.org/3/reference/datamodel.html#special-method-lookup

        """
        return cls._supress_missing


    @classmethod
    def upper(cls):
        """A reusable, but non-reentrant, context manager for temporary converting all keys passed to :meth:`__delitem__`, :meth:`__setitem__` and :meth:`__getitem__` to upper case.

        Example:

        .. code:: python

            >>> s = Settings()

            >>> s.FOO = 'bar'  # __setitem__
            >>> print(s.foo)
            foo:    bar

            >>> print(s.fOO)  # __getitem__
            bar

            >>> del s.FOO
            >>> print('FOO' in s)  # __delitem__
            False

        """
        return cls._upper


    @classmethod
    def lower(cls):
        """A reusable, but non-reentrant, context manager for temporary converting all keys passed to :meth:`__delitem__`, :meth:`__setitem__` and :meth:`__getitem__` to lower case.

        Example:

        .. code:: python

            >>> s = Settings()

            >>> s.foo = 'bar'  # __setitem__
            >>> print(s)
            FOO:    bar

            >>> print(s.FoO)  # __getitem__
            bar

            >>> del s.foo  # __delitem__
            >>> print('foo' in s)
            False

        """
        return cls._lower

    def get_nested(self, key_tuple, supress_missing=False):
        """Retrieve a nested value by, recursively, iterating through this instance using the keys in *key_tuple*.

        The :meth:`.Settings.__getitem__` method is called recursively on this instance until all keys in key_tuple are exhausted.

        Setting *supress_missing* to ``True`` will internally open the :meth:`.Settings.supress_missing` context manager, thus raising a :exc:`KeyError` if a key in *key_tuple* is absent from this instance.

        .. code:: python

            >>> s = Settings()
            >>> s.a.b.c = True
            >>> value = s.get_nested(('a', 'b', 'c'))
            >>> print(value)
            True
        """
        s = self
        with contextlib.suppress() if not supress_missing else s.supress_missing():
            for k in key_tuple:
                s = s[k]
        return s



    def set_nested(self, key_tuple, value, supress_missing=False):
        """Set a nested value by, recursively, iterating through this instance using the keys in *key_tuple*.

        The :meth:`.Settings.__getitem__` method is called recursively on this instance, followed by :meth:`.Settings.__setitem__`, until all keys in key_tuple are exhausted.


        Setting *supress_missing* to ``True`` will internally open the :meth:`.Settings.supress_missing` context manager, thus raising a :exc:`KeyError` if a key in *key_tuple* is absent from this instance.

        .. code:: python

            >>> s = Settings()
            >>> s.set_nested(('a', 'b', 'c'), True)
            >>> print(s)
            a:
              b:
                c: 	True
        """
        s = self
        *initial_keys, final_key = key_tuple
        with contextlib.suppress() if not supress_missing else s.supress_missing():
            for k in initial_keys:
                s = s[k]
        s[final_key] = value



    def flatten(self, flatten_list=True):
        """Return a flattened copy of this instance.

        New keys are constructed by concatenating the (nested) keys of this instance into tuples.

        Opposite of the :meth:`.Settings.unflatten` method.

        If *flatten_list* is ``True``, all nested lists will be flattened as well. Dictionary keys are replaced with list indices in such case.

        .. code-block:: python

            >>> s = Settings()
            >>> s.a.b.c = True
            >>> print(s)
            a:
              b:
                c: 	True

            >>> s_flat = s.flatten()
            >>> print(s_flat)
            ('a', 'b', 'c'): 	True
        """
        cls = type(self)
        if flatten_list:
            nested_type = (cls, list)
            iter_type = lambda x: x.items() if isinstance(x, cls) else enumerate(x)
        else:
            nested_type = cls
            iter_type = cls.items

        def _concatenate(key_ret, sequence):
            # Switch from Settings.items() to enumerate() if a list is encountered
            for k, v in iter_type(sequence):
                k = key_ret + (k, )
                if isinstance(v, nested_type) and v:  # Empty lists or Settings instances will return ``False``
                    _concatenate(k, v)
                else:
                    ret[k] = v

        # Changes keys into tuples
        ret = cls()
        _concatenate((), self)
        return ret



    def unflatten(self, unflatten_list=True):
        """Return a nested copy of this instance.

        New keys are constructed by expanding the keys of this instance (*e.g.* tuples) into new nested |Settings| instances.

        If *unflatten_list* is ``True``, integers will be interpretted as list indices and are used for creating nested lists.

        Opposite of the :meth:`.Settings.flatten` method.

        .. code-block:: python

            >>> s = Settings()
            >>> s[('a', 'b', 'c')] = True
            >>> print(s)
            ('a', 'b', 'c'): 	True

            >>> s_nested = s.unflatten()
            >>> print(s_nested)
            a:
              b:
                c: 	True
        """
        cls = type(self)
        ret = cls()
        for key, value in self.items():
            s = ret
            for k1, k2 in zip(key[:-1], key[1:]):
                if not unflatten_list:
                    s = s[k1]
                    continue

                if isinstance(k2, int) and not isinstance(s[k1], list):
                    s[k1] = []
                if isinstance(k1, int):  # Apply padding to s
                    s += [cls()] * (k1 - len(s) + 1)
                s = s[k1]
            s[key[-1]] = value

        return ret


    #=======================================================================


    def __iter__(self):
        """Iteration through keys follows lexicographical order. All keys are sorted as if they were strings."""
        return iter(sorted(self.keys(), key=str))


    def __missing__(self, name):
        """When requested key is not present, add it with an empty |Settings| instance as a value.

        This method is essential for automatic insertions in deeper levels. Without it things like::

            >>> s = Settings()
            >>> s.a.b.c = 12

        will not work.

        The behaviour of this method can be supressed by initializing the :class:`.Settings.supress_missing` context manager.
        """
        cls = type(self)
        self[name] = cls()
        return self[name]


    def __contains__(self, name):
        """Like regular ``__contains`__``, but if the key is an "ig" string, ignore the case."""
        if isinstance(name, ig):
            name = self.find_case(name)
        return super().__contains__(name)


    def __getitem__(self, name):
        """Like regular ``__getitem__``, but if the key is an "ig" string, ignore the case."""
        if isinstance(name, ig):
            name = self.find_case(name)
        return super().__getitem__(name)


    def __setitem__(self, name, value):
        """Like regular ``__setitem__``, but if the value is a dict, convert it to |Settings|."""
        cls = type(self)
        if isinstance(name, ig):
            name = self.find_case(name)
        if isinstance(value, dict):
            value = cls(value)
        super().__setitem__(name, value)


    def __delitem__(self, name):
        """Like regular ``__detitem__``, but if the key is an "ig" string, ignore the case."""
        if isinstance(name, ig):
            name = self.find_case(name)
        return super().__delitem__(name)


    def __getattr__(self, name):
        """If name is not a magic method, redirect it to ``__getitem__``."""
        if (name.startswith('__') and name.endswith('__')):
            return super().__getattribute__(name)
        return self[name]


    def __setattr__(self, name, value):
        """If name is not a magic method, redirect it to ``__setitem__``."""
        if name.startswith('__') and name.endswith('__'):
            super().__setattr__(name, value)
        else:
            self[name] = value


    def __delattr__(self, name):
        """If name is not a magic method, redirect it to ``__delitem__``."""
        if name.startswith('__') and name.endswith('__'):
            super().__delattr__(name)
        else:
            del self[name]


    def _str(self, indent):
        """Print contents with *indent* spaces of indentation. Recursively used for printing nested |Settings| instances with proper indentation."""
        cls = type(self)
        ret = ''
        for key, value in self.items():
            ret += ' '*indent + str(key) + ': \t'
            if isinstance(value, cls):
                ret += '\n' + value._str(indent+len(str(key))+1)
            else:  # Apply consistent indentation at every '\n' character
                indent_str = ' ' * (2 + indent + len(str(key))) + '\t'
                ret += textwrap.indent(str(value), indent_str)[len(indent_str):] + '\n'
        return ret


    def __str__(self):
        return self._str(0)

    __repr__ = __str__
    __iadd__ = soft_update
    __add__ = merge
    __copy__ = copy


class ig(str):
    """Special string that makes |Settings| work case-insensitive. Behaves exactly like the built-in `str` type. Usage: ``s = ig('abcdef')``."""
    pass
