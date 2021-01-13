#!/usr/bin/env python

"""
Class representing an XYZ MD file with changing numbers of atoms
"""

import numpy
from ..tools.periodic_table import PT
from ..mol.molecule import Molecule
from ..mol.atom import Atom
from ..core.settings import Settings
from .xyzfile import XYZTrajectoryFile
from .rkffile import bohr_to_angstrom

__all__ = ['XYZHistoryFile']

class XYZHistoryFile (XYZTrajectoryFile) :

        def __init__ (self, filename, mode='r', fileobject=None, ntap=None) :
                """
                Initializes the XYZHistoryFile object
                """
                XYZTrajectoryFile.__init__(self,filename,mode,fileobject,ntap)

                self.input_elements = self.elements[:]

        def _is_endoffile (self) :
                """
                If the end of file is reached, return coords and cell as None
                """
                end = False
                line = self.file_object.readline()
                if len(line) == 0 :
                        end = True
                        return end
                nats = int(line.split()[0])
                for i in range(nats+1) :
                        line = self.file_object.readline()
                        if len(line) == 0 :
                                end = True
                                break
                return end

        def read_coordinates (self, molecule) :
                """
                Read the coordinates at current step
                """
                # Find the number of atoms
                line = self.file_object.readline()
                if len(line) == 0 :
                        return None, None           # End of file is reached
                nats = int(line.split()[0])
                line = self.file_object.readline()

                # Read coordinates and elements
                coords = []
                elements = []
                for i in range(nats) :
                        line = self.file_object.readline() 
                        words = line.split()
                        coords.append([float(w) for w in words[1:4]])
                        elements.append(words[0])

                # If the elements changed, update the molecule
                if elements != self.elements :
                        self.elements = elements
                        self.coords = numpy.array(coords)
                        # Rebuild the molecule (bonds will disappear for now)
                        if isinstance(molecule,Molecule) :
                                for at in reversed(molecule.atoms) :
                                        molecule.delete_atom(at)
                                molecule.properties = Settings()
                                for el in elements :
                                        atom = Atom(PT.get_atomic_number(el))
                                        molecule.add_atom(atom)
                else :
                        self.coords[:] = coords

                # Assign the data to the molecule object
                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords,None,molecule,bonds=None)

        def write_next (self,coords=None,molecule=None,elements=None,cell=[0.,0.,0.],energy=0.,step=None,conect=None) :
                """             
                Write the coordinates to the next available spot in the file
                """
                if isinstance(molecule,Molecule) :
                        coords, cell, elements = self._read_plamsmol(molecule)[:3]
                self.elements = elements
                cell = self._convert_cell(cell)
                        
                self.write_moldata(coords, cell, energy, step)
                
                self.position += 1
