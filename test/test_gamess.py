
from plams import (GamessJob, Settings, finish, init)

import fnmatch
import os
import shutil


data1 = """Methylene...1-A-1 state...RHF/STO-2G
Cnv  2\n
C
H  1 rCH
H  1 rCH  2 aHCH\n
rCH=1.09
aHCH=110.0
"""


def wrapper_test():
    try:
        init()
        fun_methylene()
    finally:
        finish()
        remove_plams_dir()


def fun_methylene():
    """
    Run a RHF calculation of CH2
    """

    s = Settings()

    # data
    s.input.data = data1
    # basis
    s.input.basis.gbasis = 'sto'
    s.input.basis.ngauss = 2
    # guess
    s.input.guess.guess = 'huckel'
    # statpt
    s.input.statpt.opttol = 1e-5
    # system
    s.input.system.timlim = 1
    # control
    s.input.contrl.scftyp = 'rhf'
    s.input.contrl.runtyp = 'optimize'
    s.input.contrl.coord = 'zmt'
    s.input.contrl.nzvar = 0

    j = GamessJob(settings=s, name='methylene')
    j.run()

    r = j.results.grep_output("TOTAL ENERGY")[-1]
    final_energy = float(r.split()[-1])
    print("Total Energy: ", final_energy)

    assert abs(final_energy + 37.2380397668) < 1e-8


def remove_plams_dir(path='.'):
    """Remove Plams folder"""
    fs = os.listdir(path)
    dir_path = fnmatch.filter(fs, 'plams.*')[0]
    shutil.rmtree(dir_path)
