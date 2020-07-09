Settings
-------------

.. currentmodule:: scm.plams.core.settings

The |Settings| class provides a general purpose data container for various kinds of information that need to be stored and processed by PLAMS environment.
Other PLAMS objects (like for example |Job|, |JobManager| or |GridRunner|) have their own |Settings| instances that store data defining and adjusting their behavior.
The global scope |Settings| instance (``config``) is used for global settings.

It should be stressed here that there are no different types of |Settings| in the sense that there are no special subclasses of |Settings| for job settings, global settings etc.
Everything is stored in the same type of object and the role of a particular |Settings| instance is determined only by its content.



Tree-like structure
~~~~~~~~~~~~~~~~~~~~~~~~~

The |Settings| class is based on the regular Python dictionary (built-in class :class:`dict`, tutorial can be found :ref:`here<tut-dictionaries>`) and in many aspects works just like it::

    >>> s = Settings()
    >>> s['abc'] = 283
    >>> s[147147] = 'some string'
    >>> print(s['abc'])
    283
    >>> del s[147147]

The main difference is that data in |Settings| can be stored in multilevel fashion, whereas an ordinary dictionary is just a flat structure of key-value pairs.
That means a sequence of keys can be used to store a value.
In the example below ``s['a']`` is itself a |Settings| instance with two key-value pairs inside::

    >>> s = Settings()
    >>> s['a']['b'] = 'AB'
    >>> s['a']['c'] = 'AC'
    >>> s['x']['y'] = 10
    >>> s['x']['z'] = 13
    >>> s['x']['foo'][123] = 'even deeper'
    >>> s['x']['foo']['bar'] = 183
    >>> print(s)
    a:
      b:    AB
      c:    AC
    x:
      foo:
          123:  even deeper
          bar:  183
      y:    10
      z:    13
    >>> print(s['x'])
    foo:
        123:    even deeper
        bar:    183
    y:  10
    z:  13

So for each key the value can be either a "proper value" (string, number, list etc.) or another |Settings| instance that creates a new level in the data hierarchy.
That way similar information can be arranged in subgroups that can be copied, moved and updated together.
It is convenient to think of a |Settings| object as a tree.
The root of the tree is the top instance (``s`` in the above example), "proper values" are stored in leaves (a leaf is a childless node) and internal nodes correspond to nested |Settings| instances (we will call them *branches*).
Tree representation of ``s`` from the example above is illustrated on the following picture:

.. image:: ../_static/set_tree.*


Tree-like structure could also be achieved with regular dictionaries, but in a rather cumbersome way::

    >>> d = dict()
    >>> d['a'] = dict()
    >>> d['a']['b'] = dict()
    >>> d['a']['b']['c'] = dict()
    >>> d['a']['b']['c']['d'] = 'ABCD'
    ===========================
    >>> s = Settings()
    >>> s['a']['b']['c']['d'] = 'ABCD'

In the last line of the above example all intermediate |Settings| instances are created and inserted automatically.
Such a behavior, however, has some downsides -- every time you request a key that is not present in a particular |Settings| instance (for example as a result of a typo), a new empty instance is created and inserted as a value of this key.
This is different from dictionaries where exception is raised in such a case::

    >>> d = dict()
    >>> d['foo'] = 'bar'
    >>> x = d['fo']
    KeyError: 'fo'
    ===========================
    >>> s = Settings()
    >>> s['foo'] = 'bar'
    >>> x = s['fo']

    >>> print(s)
    fo:            #the value here is an empty Settings instance
    foo:    bar


.. _dot-notation:

Dot notation
~~~~~~~~~~~~~~~~~~~~~~~~~

To avoid inconvenient punctuation, keys stored in |Settings| can be accessed using the dot notation in addition to the usual bracket notation.
In other words ``s.abc`` works as a shortcut for ``s['abc']``.
Both notations can be used interchangeably::

    >>> s.a.b = 'AB'
    >>> s['a'].c = 'AC'
    >>> s.x['y'] = 10
    >>> s['x']['z'] = 13
    >>> s['x'].foo[123] = 'even deeper'
    >>> s.x.foo.bar = 183
    >>> print(s)
    a:
      b:    AB
      c:    AC
    x:
      foo:
          123:  even deeper
          bar:  183
      y:    10
      z:    13

Due to the internal limitation of the Python syntax parser, keys other than single word strings cannot work with that shortcut, for example::

    >>> s.123.b.c = 12
    SyntaxError: invalid syntax
    >>> s.q we.r.t.y = 'aaa'
    SyntaxError: invalid syntax
    >>> s.5fr = True
    SyntaxError: invalid syntax

In those cases one has to use the regular bracket notation::

    >>> s[123].b.c = 12
    >>> s['q we'].r.t.y = 'aaa'
    >>> s['5fr'] = True

The dot shortcut does not work for keys which begin and end with two (or more) underscores (like ``__key__``).
This is done on purpose to ensure that Python magic methods work properly.



Case sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~

|Settings| are case-preserving but case-insensitive. That means every key is stored in its original form, but when looked up (for example, to access value, test existence or delete), any casing can be used::

    >>> s = Settings()
    >>> s.foo = 'bar'
    >>> s.System.one = 1
    >>> s.system.two = 2
    >>> print(s.FOO)
    bar
    >>> 'Foo' in s
    True
    >>> print(s)
    foo:    bar
    System:
       one:     1
       two:     2
    >>> 'oNe' in s.SYSTEM
    True



Global settings
~~~~~~~~~~~~~~~~~~~~~~~~~

Global settings are stored in a public |Settings| instance named ``config``.
They contain variables adjusting general behavior of PLAMS as well as default settings for various objects (jobs, job manager etc.)
The ``config`` instance is created during initialization of PLAMS environment (see |init|) and populated by executing ``plams_defaults`` file.
It is visible in the main PLAMS namespace so every time you wish to adjust some settings you can simply type in your script, for example::

    config.job.pickle = False
    config.sleepstep = 10

These changes are going to affect only the script they are called from.
If you wish to permanently change some setting for all PLAMS executions, you can do it by editing ``plams_defaults``, which is located in the root folder of the package (``$AMSHOME/scripting/scm/plams``).

.. note::

    You can create multiple "profiles" of PLAMS behavior by creating multiple different copies of ``plams_defaults`` (also with different filenames).
    If the environmental variable ``$PLAMSDEFAULTS`` is present and its value points to an existing file, this file is used instead of ``plams_defaults`` from the root folder.



API
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: Settings
    :exclude-members: __weakref__, __copy__, __add__, __iadd__

.. note::

    Methods :meth:`~Settings.update` and :meth:`~Settings.soft_update` are complementary.
    Given two |Settings| instances ``A`` and ``B``, the command ``A.update(B)`` would result in ``A`` being exactly the same as ``B`` would be after ``B.soft_update(A)``.
