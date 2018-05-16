.. _public-functions:

Public functions
-------------------------

.. currentmodule:: scm.plams.core.functions

This chapter gathers information about publicly visible functions that can be used in PLAMS scripts.

.. autofunction:: init
.. autofunction:: finish
.. autofunction:: load
.. autofunction:: load_all
.. autofunction:: scm.plams.tools.load.load_manual
.. autofunction:: scm.plams.tools.load.load_all_manual
.. autofunction:: scm.plams.core.basemol.read_molecules

.. _logging:

Logging
~~~~~~~~~~~~~~~~~~~~~~~~~

PLAMS features a simple logging mechanism. All important actions happening in functions and methods register their activity using log messages. These massages can be printed to standard output and/or saved to ``.log`` file located in the main working folder.

Every log message has its "verbosity" defined as an integer number: the higher the number, the more detailed and descriptive the message is. In other words, it is a measure of importance of the message. Important events (like "job started", "job finished", "something went wrong") should have low verbosity, whereas less crucial ones (for example "pickling of job X successful") a bit higher. The purpose of that is to allow user to choose how verbose the whole logfile is. Each logfile (either file or stdout) has an integer number associated with it defining which messages are printed to this logfile (for example, if this number is 3, all messages with verbosity 3 or less are printed). That way picking a smaller number results in logfile being short and containing only the most relevant information while larger numbers produce longer and more detailed logfiles.

The behavior of the logging mechanism is adjusted by ``config.log.`` settings branch with the following keys:

*   ``file`` (integer) -- verbosity of logfile printed to the ``.log`` file.
*   ``stdout`` (integer) -- verbosity of logfile printed to the standard output.
*   ``time`` (boolean) -- print time of each log event.
*   ``date`` (boolean) -- print date of each log event.

Log messages used within the PLAMS code use four different levels of verbosity:

*   **1**: important
*   **3**: normal
*   **5**: verbose
*   **7**: debug

Even levels are left empty for user's convenience. For example, if you find level 5 too verbose and still want to be able to switch on and off log messages of your own code, you can log them with verbosity 4.

.. note::

    Your own code can (and should) contain some |log| calls. They are very important for debugging purposes.

.. autofunction:: log



.. _binding-decorators:

Binding decorators
~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes one wants to expand functionality of a class by adding a new method or modifying an existing one. It can be done in a few different ways:

*   One can go directly to the source code defining the class and modify it there before running a script. Such a change is global -- it affects all the future scripts, so in most cases it is not a good thing (for defining |prerun| for example).
*   Creating a subclass with new/modified methods definitions is usually the best solution. It can be done directly in your script before the work is done or in a separate dedicated file executed before the actual script (see :ref:`master-script`). Newly defined class can be then used instead of the old one and changes will be reflected. However, this solution fails in some rare cases when a method needs to differ for different instances or when it needs to be changed during the runtime of the script.
*   PLAMS binding decorators (|add_to_class| and |add_to_instance|) can be used.

Binding decorators allow to bind methods to existing classes or even directly to particular instances without having to define a subclass. Such changes are visible only inside the script in which they are used.

To fully understand how binding decorators work let us take a look at how Python handles method calling. Assume we have an instance of a class (let's say ``myres`` is an instance of |DFTBResults|) and there is a method call in our script (let it be ``myres.somemethod(arguments)``). Python first looks for ``somemethod`` amongst attributes of ``myres``. If it is not there (which is usually the case, since methods are defined in classes), attributes of |DFTBResults| class are checked. If ``somemethod`` is still not there, parent classes are checked in the order of inheritance (in our case, first |SCMResults|, then |Results|). That implies two important things:

*   |add_to_instance| affects only one particular instance, but is "stronger" than |add_to_class| -- method added to instance always takes precedence before the same method added to (or just defined in) a class
*   changes done with |add_to_class| affect all instances of that particular class, including even those created before |add_to_class| was used.

The usage of binding decorators is straightforward. You simply define a regular function somewhere inside your script and decorate it with one of the decorators (see below). The function needs to be a valid method code, so it should have ``self`` as the first argument and use it to reference the class instance.


.. autofunction:: add_to_class
.. autofunction:: add_to_instance

.. technical::

    Each of the above decorators is in fact a decorator factory that, given an object (class or instance), produces a decorator that binds function as a method of that object. Both decorators are adding instance methods only, they cannot be used for static or class methods.
