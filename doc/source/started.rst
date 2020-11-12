Getting started
=========================

This section contains general information about installing and running PLAMS.



Library contents
-------------------------

PLAMS comes in a form of Python 3 package (earlier versions were also compatible with Python 2, but as Python 3 becomes more and more popular, we decided to drop Python 2 compatibility).
The root folder of the package contains the following:

*   ``core``: core subpackage containing all the essential modules defining the skeleton of the library
*   ``interfaces``: subpackage with interfaces to external binaries
*   ``tools``: subpackage with small utilities like unit converter, periodic table, file readers etc.
*   ``recipes``: subpackage with examples of simple job types built on top basic PLAMS elements
*   ``doc``: Sphinx source of this documentation
*   ``scripts``: folder with executable scripts
*   ``plams_defaults``: a separate file defining all the adjustable options accessible via ``config`` object (see below for details)


An imporant part of PLAMS worth mentioning here is the executable script used for running or restarting your workflows.
It is called simply ``plams`` and it's located in the ``scripts`` folder.
We will refer to it as the *launch script* or simply *launcher*.
Further in this section you can find the dedicated chapter explaining the usage of the launch script.



Installing PLAMS
-------------------------

You can install PLAMS on your computer using one of the following ways:

1.  If you are using Amsterdam Modeling Suite, PLAMS is shipped as a part of ``scm`` Python package (``$AMSHOME/scripting/scm/plams``) and configured to work with a built-in Python coming with AMSuite (you can access it with ``amspython`` command).
    The launch script is added to ``$AMSBIN``, so it should be directly visible from your command line (as long as ``$AMSBIN`` is in your ``$PATH``).

2.  The latest PLAMS stable release can be installed directly from PyPi by typing ``pip install plams`` in your command line.
    The launch scipt will be installed along other global system executables (platform dependent) and should be visible from your command line.

3.  Any current or historic version can be downloaded or cloned from PLAMS `GitHub page <https://github.com/SCM-NV/PLAMS>`_.
    The ``release`` branch points to the latests stable release, while the ``master`` branch is the most recent development snapshot.

4.  You can combine methods 2 and 3 and fetch PLAMS from GitHub using ``pip``: ``pip install git+https://github.com/SCM-NV/PLAMS.git@master`` (make sure to have Git installed and to choose the proper branch)

PLAMS requires the following Python packages as dependencies:

*   `numpy <http://www.numpy.org>`_
*   `dill <https://pypi.python.org/pypi/dill>`_ (enhanced pickling)
*   `ase <https://wiki.fysik.dtu.dk/ase>`_ (optional dependency)
*   `rdkit <https://pypi.org/project/rdkit>`_ (optional dependency)

If you are using Amsterdam Modeling Suite, all the above packages are already included in our Python stack.
When you install PLAMS using ``pip``, the required packages (numpy and dill)will be installed automatically.
For optional dependencies, or in any other case you can install them with ``pip install [package name]``.



Updating PLAMS
-------------------------

PLAMS is still a relatively new project and new developments based on users' feedback are added frequently.
Because of that some users might prefer to use the recent development version (``master`` branch on GitHub) rather than the usual stable release (``release`` branch).

If you're running PLAMS as a part of Amsterdam Modeling Suite you can update the library to the recent development snapshot simply by downloading a ZIP archive from PLAMS `GitHub page <https://github.com/SCM-NV/PLAMS>`_ and replacing contents of ``$AMSHOME/scripting/scm/plams`` with the contents of that archive.

If you obtained PLAMS using Git or ``pip`` please use the usual Git or ``pip`` methods to pull a new version (``git pull`` or ``pip install --upgrade plams``)



Running PLAMS
-------------------------

Inside your Python interpreter or in Python scripts PLAMS is visible as a subpackage of the ``scm`` package, so you can import it with one of the following commands::

    import scm.plams
    from scm import plams
    from scm.plams import *

.. note::

    Usually in Python ``import *`` is considered a bad practice and discouraged.
    However, PLAMS internally takes care of the namespace tidiness and imports only necessary things with ``import *``.
    Importing with ``import *`` allows you to use identifiers like ``Molecule`` or ``AMSJob`` instead of ``scm.plams.Molecule`` or ``scm.plams.AMSJob`` which makes your scripts shorter and more readable.
    Throughout this documentation it is assumed that ``import *`` is used so identifiers are not prefixed with ``scm.plams.`` in any example.

.. technical::

    PLAMS namespace is automatically built when importing the package.
    All Python modules present in any of four main subpackages (``core``, ``tools``, ``interfaces``, ``recipes``) are processed and each module's ``__all__`` attributes are added to the main namespace.
    That scheme works in plug'n'play manner -- if you add a new ``.py`` file with properly defined ``__all__`` attribute in one of the four abovementioned folders, names defined in that ``__all__`` attribute will be added to the main namespace.


A PLAMS script is in fact a general Python script that makes use of classes and functions defined in the PLAMS library.
To work properly such a script has to follow two simple restrictions.
At the very beginning of the script one must initialize PLAMS environment by calling public function |init|.
Without this initialization almost every PLAMS function or class call results in a crash.
Similarly, at the end of the script public function |finish| needs to be called to properly clean the main working folder and ensure proper closure of parallel scripts.
You can find more detailed information about these two functions in :ref:`public-functions` section.

To sum up, a proper PLAMS script needs to look like this::

    from scm.plams import *
    init()
    # =========
    # actual script here
    # ...
    # =========
    finish()

and it should be executed from the command line with ``python [filename]`` (``amspython [filename]`` in case of AMS Python stack).
Keeping these restrictions in mind can be a bit inconvenient, so PLAMS comes with the launcher that takes care of the proper initialization and cleaning.
See |master-script| for details.

Of course PLAMS can be also run interactively.
After starting your favorite Python interpreter you need to manually import and initialize the environment with ``from scm.plams import *`` and |init|.
Then you can interactively run any Python command relying on PLAMS.
If you run any jobs in the interactive mode make sure to use |finish| before closing the interpreter to ensure that all the jobs are gently finished and the main working folder is cleaned.



.. _plams-defaults:

Defaults file
-------------------------

The defaults file is called ``plams_defaults`` and it is located in the root folder of the package.
If you installed PLAMS using ``pip``, the defaults file could be a bit difficult to find (usually somewhere in ``site-packages`` subfolder of your Python).
If you can't find it, just get a fresh copy from `GitHub <https://github.com/SCM-NV/PLAMS/blob/master/plams_defaults>`_, put it somewhere on your disk and set ``$PLAMSDEFAULTS`` environmental variable pointing to it.
See also |init| to find how PLAMS looks for the defaults file.

The defaults file contains a list of commands that adjust various aspects of PLAMS behavior.
The file is self-explanatory: each command is preceded with a comment explaining what it does.
We **strongly recommend** to have a quick glance at that file.
It gives an overview of what and how can be tweaked (it's not long, we promise).

If you wish to globally change some setting you can do it by modifying the defaults file.
Changes you make there are going to affect all future PLAMS runs.
To tweak a particular setting just for a single script, copy a corresponding line from the defaults file and place it at the top of your script.
For example::

    config.log.stdout = 1
    config.job.pickle = False
    config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)



.. _master-script:

The launch script
-------------------------

The launch script is an executable file called simply ``plams`` located in the ``scripts`` folder.
If your ``$PATH`` variable is configured properly, you can type in your command line ``plams -h`` or ``plams --help`` for a short help message.

The launch script provides a convenient way of executing PLAMS scripts and takes care of important things mentioned earlier in this chapter: properly importing and initializing PLAMS and cleaning after all the work is done.
Thanks to that your actual script does not need to contain import, init or finish commands.

Without the launcher::

    from scm.plams import *
    init()
    # =========
    # actual script here
    # ...
    # =========
    finish()

executed with ``python [filename]`` (or ``amspython [filename]``).

With the launcher::

    # =========
    # actual script here
    # ...
    # =========

executed with ``plams [filename]``.

Besides that, the launch script offers several command line arguments allowing you to tune the behavior of your script without a need to edit the script itself.


Working folder location
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script allows you to pick custom name and location for the main working folder.
The main working folder is an initially empty folder that is created on |init|.
All files produced by PLAMS and other programs executed by it are saved in the main working folder (usually in some of its subfolders).
Each separate run of PLAMS has its separate main working folder.

By default the main working folder is located in the directory where your script was executed and is called ``plams_workdir`` (``plams_workdir.002`` if ``plams_workdir`` already existed).
You can change that by supplying ``-p`` and ``-f`` (or ``--path`` and ``--folder``) arguments to the launcher to choose the location and the name of the main working folder.
For example the command::

    plams -p /home/user/science -f polymers myscript.plms

will use ``/home/user/science/polymers`` as the main working folder regardless where this command was executed.

.. note::

    Each PLAMS run creates a fresh, empty directory for its main working folder.
    If you try to use an existing folder (or don't pick any and ``plams_workdir`` already exists in the current directory), a unique folder is going to be created anyway, by appending ``.002`` (or ``.003``, ``.004`` and so on) to the name of your folder.


Passing variables
~~~~~~~~~~~~~~~~~~~~~~~~~

When using the launcher you can pass variables to your script directly from the command line.
This can be done with ``-v`` (or ``--var``) parameter that follows the syntax ``-v variable=value`` (mind the lack of spaces around equal sign, it is a must).
For a script executed that way, there is an additional global string variable with the name ``variable`` and the value ``'value'`` visible from within the script.
For example if the script in file ``script1.plms`` looks like this::

    print('Chosen basis: ' + basis)
    print('Number of points: ' + n)
    print(type(n))
    # do something depending on n and basis

and you execute it with::

    plams -v n=10 -v basis=DZP script1.plms

the standard output will be:

.. code-block:: none

    Chosen basis: DZP
    Number of points: 10
    str
    [rest of the output]

Three important things to keep in mind about ``-v`` parameter:

*   no spaces around equal sign,
*   each variable requires separate ``-v``,
*   the type of the variable is **always** string (like in the example above).
    If you want to pass some numerical values, make sure to convert them from strings to numbers inside your script.


Importing past jobs
~~~~~~~~~~~~~~~~~~~~~~~~~

You can instruct the launcher to load the results of some previously run jobs by supplying the path to the main working folder of a finished PLAMS run with ``-l`` (or ``--load``) parameter.
To find out why this could be useful, please see |pickling| and |RPM|.

This mechanism is equivalent to using |load_all| function at the beginning of your script.
That means executing your script with ``plams -l /some/path myscript.plms`` works just like putting ``load_all('/some/path')`` at the beginning of ``myscript.plms`` and running it with ``plams myscript.plms``.
The only difference is that, when using |load_all| inside the script, you can access each of the loaded jobs separately by using the dictionary returned by |load_all|.
This is not possible with ``-l`` parameter, but all the loaded jobs will be visible to |RPM|.

Multiple different folders can be supplied with ``-l`` parameter, but each of them requires a separate ``-l`` flag::

    plams -l /some/path -l /other/path myscript.plms


Restarting failed script
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script can be called with an additional argumentless ``-r`` parameter (or ``--restart``).
In such a case the launcher enters "restart mode".
In the restart mode the folder specified by ``-f`` (or the latest ``plams_workdir[.xxx]`` if ``-f`` is not used) is first renamed by appending ``.res`` to folder's original name (let's call it ``foldername``).
Successful jobs from ``foldername.res`` are loaded at the beginning of the current run, which is executed in a new, empty main working folder called ``foldername``.
Whenever the new run encounters a job identical to a successful job present in ``foldername.res``, the new job execution is skipped and the whole job folder is linked (hardlinked) from ``foldername.res`` to ``foldername``.
That way the restart run will not redo any work present in old ``foldername``, but rather back it up to ``foldername.res`` and restart from the point when the old run was terminated.
For example, after::

    $ plams -f stuff myscript.plms
    [17:28:40] PLAMS working folder: /home/user/stuff
    # [some successful work]
    [17:56:22] Execution interrupted by the following exception:
    # [exception details]

you can edit ``myscript.plms``, remove the cause of crash and restart your script with::

    $ plams -f stuff -r myscript.plms
    RESTART: Moving stuff to stuff.res and restarting from it
    [18:03:34] PLAMS working folder: /home/user/stuff

(the above command needs to be executed in ``/home/user``.
Otherwise, you need to add ``-p /home/user`` to tell the master script where to look for ``stuff``).
The same example with the default folder name::

    $ plams myscript.plms
    [17:28:40] PLAMS working folder: /home/user/plams_workdir
    # [some successful work]
    [17:56:22] Execution interrupted by the following exception:
    # [exception details]

    [...debug the script...]

    $ plams -r myscript.plms
    RESTART: Moving plams_workdir to plams_workdir.res and restarting from it
    [18:03:34] PLAMS working folder: /home/user/plams_workdir

For more detailed explanation of the restart mechanism, please see |RPM|, |pickling| and |restarting|.


Multiple input scripts
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script can be called with more than one positional argument, like for example::

    plams script1.plms script2.plms script3.plms

All files supplied that way are concatenated into one script and then executed (that means things declared in script1 are visible in script2 and script3).
Using this feature for completely unrelated scripts is probably not a good idea, but it can be useful, for example, when first files contain just definitions of your own functions, derived classes, settings tweaks etc. that are then used in the last file::

    plams config/debug_run.plms settings/adf/adf_fde.plms actual_script.plms

That way you can build your own library of reusable code snippets for tasks that are most frequently occurring in your daily work, customize PLAMS according to your personal preferences and make your working environment truly modular.

.. note::

    The ``.plms`` file extension for PLAMS scripts is just a convention.
    Scripts can be any text files.

