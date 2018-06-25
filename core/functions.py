import builtins
import glob
import os
import shutil
import sys
import threading
import time
import types

from os.path import join as opj
from os.path import isfile, isdir, expandvars, dirname

from .errors import PlamsError
from .settings import Settings

__all__ = ['init', 'finish', 'log', 'load', 'load_all', 'add_to_class', 'add_to_instance']


#===========================================================================


def init(path=None, folder=None):
    """Initialize PLAMS environment. Create global ``config`` and default |JobManager|.

    An empty |Settings| instance is created and added to :mod:`public<builtins>` namespace as ``config``. Then it is populated with default settings by executing ``plams_defaults``. The following locations are used to search for the defaults file, in order of precedence:

    *   If ``$PLAMSDEFAULTS`` variable is in your environment and it points to a file, this file is used (executed as Python script).
    *   If ``$PLAMSHOME`` variable is in your environment and ``$PLAMSHOME/src/scm/plams/plams_defaults`` exists, it is used.
    *   If ``$ADFHOME`` variable is in your environment and ``$ADFHOME/scripting/plams/src/scm/plams/plams_defaults`` exists, it is used.
    *   Otherwise, the path ``../../plams_defaults`` relative to the current file (``functions.py``) is checked. If defaults file is not found there, an exception is raised.

    Next, a |JobManager| instance is created as ``config.jm`` using *path* and *folder* to determine the main working directory. Settings used by this instance are directly linked from ``config.jobmanager``. If *path* is not supplied, the current directory is used. If *folder* is not supplied, the string ``plams.`` followed by PID of the current process is used.

    .. warning::
      This function **must** be called before any other PLAMS command can be executed. Trying to do anything without it results in a crash. See also :ref:`master-script`.
    """

    builtins.config = Settings()

    if 'PLAMSDEFAULTS' in os.environ and isfile(expandvars('$PLAMSDEFAULTS')):
        defaults = expandvars('$PLAMSDEFAULTS')
    elif 'PLAMSHOME' in os.environ and isfile(opj(expandvars('$PLAMSHOME'), 'src', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$PLAMSHOME'), 'src', 'scm', 'plams', 'plams_defaults')
    elif 'ADFHOME' in os.environ and isfile(opj(expandvars('$ADFHOME'), 'scripting', 'plams', 'src', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$ADFHOME'), 'scripting', 'plams', 'src', 'scm', 'plams', 'plams_defaults')
    else:
        defaults = opj(dirname(dirname(__file__)), 'plams_defaults')
        if not isfile(defaults):
            raise PlamsError('plams_defaults not found, please set PLAMSDEFAULTS or PLAMSHOME in your environment')
    exec(compile(open(defaults).read(), defaults, 'exec'))

    from .jobmanager import JobManager
    config.jm = JobManager(config.jobmanager, path, folder)

    log('Running PLAMS located in {}'.format(dirname(dirname(__file__))) ,5)
    log('Using Python {}.{}.{} located in {}'.format(*sys.version_info[:3], sys.executable), 5)
    log('PLAMS defaults were loaded from {}'.format(defaults) ,5)

    log('PLAMS environment initialized', 5)
    log('PLAMS working folder: {}'.format(config.jm.workdir), 1)

    try:
        import dill
    except ImportError:
        log('WARNING: importing dill failed. Falling back to default pickle module. Expect problems with pickling', 1)


#===========================================================================


def finish(otherJM=None):
    """Wait for all threads to finish and clean the environment.

    This function must be called at the end of your script for :ref:`cleaning` to take place. See :ref:`master-script` for details.

    If for some reason you use other job managers than the default one, they need to passed as *otherJM* list.
    """
    for thread in threading.enumerate():
        if thread.name == 'plamsthread':
            thread.join()

    config.jm._clean()
    if otherJM:
        for jm in otherJM:
            jm._clean()
    log('PLAMS environment cleaned up successfully', 5)
    log('PLAMS run finished. Goodbye', 3)

    if config.erase_workdir is True:
        shutil.rmtree(config.jm.workdir)


#===========================================================================


def load(filename):
    """Load previously saved job from ``.dill`` file. This is just a shortcut for |load_job| method of the default |JobManager| ``config.jm``."""
    return config.jm.load_job(filename)


#===========================================================================


def load_all(path, jobmanager=None):
    """Load all jobs from *path*.

    This function works as a multiple execution of |load_job|. It searches for ``.dill`` files inside the directory given by *path*, yet not directly in it, but one level deeper. In other words, all files matching ``path/*/*.dill`` are used. That way a path to the main working folder of a previously run script can be used to import all the jobs run by that script.

    In case of partially failed |MultiJob| instances (some children jobs finished successfully, but not all) the function will search for ``.dill`` files in children folders. That means, if ``path/[foldername]/`` contains some subfolders (for children jobs) but does not contail a ``.dill`` file (the |MultiJob| was not fully successful), it will look into these subfolders. This behavior is recursive up to arbitrary folder tree depth.

    The purpose of this function is to provide quick and easy way of restarting a script that previously failed. Loading all successful jobs from the previous run prevents double work and allows the script to proceed directly to the place where it failed.

    Jobs are loaded using default job manager stored in ``config.jm``. If you wish to use a different one you can pass it as *jobmanager* argument of this function.

    Returned value is a dictionary containing all loaded jobs as values and absolute paths to ``.dill`` files as keys.
    """
    jm = jobmanager or config.jm
    loaded_jobs = {}
    for foldername in filter(lambda x: isdir(opj(path,x)), os.listdir(path)):
        maybedill = opj(path,foldername,foldername+'.dill')
        if isfile(maybedill):
            job = jm.load_job(maybedill)
            if job:
                loaded_jobs[os.path.abspath(maybedill)] = job
        else:
            loaded_jobs.update(load_all(path=opj(path,foldername), jobmanager=jm))
    return loaded_jobs


#===========================================================================


_stdlock = threading.Lock()
_filelock = threading.Lock()

def log(message, level=0):
    """Log *message* with verbosity *level*.

    Logs are printed independently to both text file and standard output. If *level* is equal or lower than verbosity (defined by ``config.log.file`` or ``config.log.stdout``) the message is printed. Date and/or time can be added based on ``config.log.date`` and ``config.log.time``. All logging activity is thread safe.
    """
    if 'config' in vars(builtins):
        if level <= config.log.file or level <= config.log.stdout:
            message = str(message)
            prefix = ''
            if config.log.date:
                prefix += '%d.%m|'
            if config.log.time:
                prefix += '%H:%M:%S'
            if prefix:
                prefix = '[' + prefix.rstrip('|') + '] '
                message = time.strftime(prefix) + message
            if level <= config.log.stdout:
                with _stdlock:
                    print(message)
            if level <= config.log.file and 'jm' in config:
                with _filelock, open(config.jm.logfile, 'a') as f:
                    f.write(message + '\n')


#===========================================================================


def add_to_class(classname):
    """Add decorated function as a method to the whole class *classname*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        @add_to_class(ADFResults)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

    After executing the above code all instances of ``ADFResults`` (even the ones created earlier) are enriched with ``get_energy`` method that can be invoked by::

        someadfresults.get_energy()

    The added method is visible from subclasses of *classname* so ``@add_to_class(Results)`` will also work in the above example.

    If *classname* is |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see :ref:`parallel`).
    """
    from .results import _restrict, _MetaResults
    def decorator(func):
        if isinstance(classname, _MetaResults):
            func = _restrict(func)
        setattr(classname, func.__name__, func)
    return decorator


#===========================================================================


def add_to_instance(instance):
    """Add decorated function as a method to one particular *instance*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        results = myjob.run()

        @add_to_instance(results)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

        results.get_energy()

    The added method is visible only for one particular instance and it overrides any methods defined on class level or added with :func:`add_to_class` decorator.

    If *instance* is an instance of |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see :ref:`parallel`).
    """
    from .results import _restrict, Results
    def decorator(func):
        if isinstance(instance, Results):
            func = _restrict(func)
        func = types.MethodType(func, instance)
        setattr(instance, func.__func__.__name__, func)
    return decorator

