import copy
import hashlib
import subprocess
import time


from .functions import config, log


__all__ = []



def smart_copy(obj, owncopy=[], without=[]):
    """Return a copy of *obj*. Attributes of *obj* listed in *without* are ignored. Attributes listed in *owncopy* are copied by calling their own ``copy()`` methods. All other attributes are copied using :func:`copy.deepcopy`."""

    ret = obj.__class__()
    for k in owncopy:
        ret.__dict__[k] = obj.__dict__[k].copy()
    for k in obj.__dict__:
        if k not in (without + owncopy):
            ret.__dict__[k] = copy.deepcopy(obj.__dict__[k])
    return ret


#===========================================================================


def sha256(string):
    """A small utility wrapper around :ref:`hashlib.sha256<hash-algorithms>`."""
    if not isinstance(string, bytes):
        string = str(string).encode()
    h = hashlib.sha256()
    h.update(string)
    return h.hexdigest()


#===========================================================================


def saferun(*args, **kwargs):
    """A wrapper around :func:`subprocess.run` repeating the call ``config.saferun.repeat`` times with ``config.saferun.delay`` interval in case of :exc:`BlockingIOError` being raised (any other exception is not caught and directly passed above). All arguments (*args* and *kwargs*) are passed directly to :func:`~subprocess.run`. If all attempts fail, the last raised :exc:`BlockingIOError` is reraised."""
    attempt = 0
    (repeat, delay) = (config.saferun.repeat, config.saferun.delay) if ('saferun' in config) else (5,1)
    while attempt <= repeat:
        try:
            return subprocess.run(*args, **kwargs)
        except BlockingIOError as e:
            attempt += 1
            log('subprocess.run({}) attempt {} failed with {}'.format(args[0], attempt, e), 5)
            last_error = e
            time.sleep(delay)
    raise last_error


