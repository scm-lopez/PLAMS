import glob
import os
import shutil

from os.path import join as opj
from os.path import isfile, isdir, abspath, split, basename

from ..core.errors import PlamsError, FileError
from ..core.settings import Settings
from ..interfaces.adfsuite.adf import ADFJob

__all__ = ['load_manual', 'load_all_manual']



def load_manual(path, jobType=ADFJob, settings=None):
    """
    Creates a PLAMS Job of ``jobType`` from an existing folder, returns a |Job| object.

    * Folder Name is the job name
    * Dependencies and all Input are not recreated in PLAMS
    * Settings are set to what you pass as ``settings``, default: empty ``Settings()``
    * ``job.run = None`` to prevent usage.
    * PLAMS naming scheme has to be followed for the job files, all job files as defined by the Job class need to be present.

    Only use this to access results from Jobs run outside PLAMS.
    """
    if settings is None:
        settings = Settings()

    if not isdir(path):
        raise FileError('Path {:} does not exist, cannot load from it.'.format(path))

    if not hasattr(jobType,'_filenames'):
        raise PlamsError('You are requesting to load from jobType {:}, but it does not provide `_filenames`'.format(str(jobType)))

    oldPath = abspath(path)
    fileExt =  jobType._filenames
    name = split(oldPath)[1]

    for key in fileExt:
        f = opj(oldPath,fileExt[key]).replace('$JN',name)
        if not isfile(f):
            raise FileError('File {:} not present, cannot load job.'.format(f))

    #disable hashing for this one
    tmpSett = settings.copy()
    tmpSett.hashing = False
    tmpSett.keep = "all"
    job = jobType(name=name, settings=tmpSett)
    #job.path should now be none, check and create folder
    if job.path is None:
        job.path = opj(config.jm.workdir, job.name)
        if os.path.exists(job.path):
            if config.jm.settings.jobfolder_exists == 'remove':
                shutil.rmtree(job.path)
            elif config.jm.settings.jobfolder_exists == 'rename':
                i = 1
                while os.path.exists(job.path + '.old' + str(i)):
                    i += 1
                newname = job.path + '.old' + str(i)
                os.rename(job.path, newname)
                log('Folder {} already present. Renaming it to {}'.format(job.path, newname), 1)
        os.mkdir(job.path)
    else:
        raise PlamsError("Job {:} has just been created and already has a path.".format(job.name))

    #avoid running this job
    job.run = None

    newPath = opj(job.path,'')
    for f in glob.glob(opj(oldPath,'*')):
        if isdir(f):
            shutil.copytree(f, opj(newPath,basename(f)))
        else:
            shutil.copy(f, newPath)


    #finalize job (includes collecting results into job)
    job._finalize()

    return job


def load_all_manual(path, ignore=None, jobType=ADFJob, settings=None):
    """
    Load jobs executed without PLAMS from existing subdirectories. 
    Calls the |load_manual| function on all subdirectories not in ``ignore``.
    Ignore can be a string or list containing the folder names, not paths. Regular Expressions are treated using re.

    All Jobs need to be of the same ``jobType``!
    """

    folders = [ f for f in glob.glob(opj(path,'*')) if isdir(f) ]

    if isinstance(ignore, list):
        ign = []
        for f in ignore:
            ign.extend(glob.glob(opj(path,f)))
    else:
        ign = glob.glob(opj(path,f))

    jobs = []

    for folder in folders:
        #skip folders in ignore
        if folder in ign:
            continue

        jobs.append(load_manual(folder, jobType=jobType, settings=settings))

    return jobs
