"""
What follows is a somewhat hacky way to set up an automatic plug'n'play import mechanism. If you experience any problems with importing, remove this and bring back good old explicit imports:

        from .core.basejob import *
        from .core.basemol import *
        ...
        from .tools.geometry import *
        from .tools.kftools import *
        ...
        from .interfaces.adfsuite import *
        from .interfaces.cp2k import *
        from .interfaces.crystal import *
        from .interfaces.dftbplus import *
"""

def __autoimport(path, folders):
    import os
    from os.path import join as opj
    is_module = lambda x: x.endswith('.py') and not x.startswith('__init__')

    ret = []
    for folder in folders:
        for dirpath, dirnames, filenames in os.walk(opj(path,folder)):
            modules = [os.path.splitext(f)[0] for f in filter(is_module, filenames)]
            relpath = os.path.relpath(dirpath, path).split(os.sep)
            for module in modules:
                imp = '.'.join(relpath + [module])
                tmp = __import__(imp, globals=globals(), fromlist=['*'], level=1)
                if hasattr(tmp, '__all__'):
                    ret += tmp.__all__
                    for name in tmp.__all__:
                        globals()[name] = vars(tmp)[name]
    return ret


__all__ = __autoimport(__path__[0], ['core', 'interfaces', 'tools', 'recipes'])
