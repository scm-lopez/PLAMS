#!/usr/bin/env python

import numpy
from scm.plams import Molecule, cell_shape
from .trajectoryfile import TrajectoryFile

__all__ = ['XYZTrajectoryFile','create_xyz_string']

class XYZTrajectoryFile (TrajectoryFile) :
        """
        Class that represents an XYZ trajectory file
        """
        def __init__ (self, filename, mode='r', fileobject=None, ntap=None) :
                """
                Creates an XYZ trajectory file object

                @param ntap :     Number of atoms (constant throughout trajectory)
                @param elements : The element names of the atoms
                """
                self.position = 0
                if filename is not None :
                        fileobject = open(filename,mode)
                self.file_object = fileobject
                if self.file_object is not None :
                        self.mode = self.file_object.mode
                
                self.ntap = 0
                if ntap is not None :
                        self.ntap = ntap
                self.firsttime = True
                self.coords = numpy.zeros((self.ntap,3))                # Only for reading purposes, to avoid creating the array each time

                # XYZ specific attributes
                self.elements = ['H']*self.ntap
                self.name = 'PlamsMol'

                # Required setup before frames can be read/written
                if self.mode == 'r' :
                        self.read_header()

        def get_elements (self) :        
                """
                Get the elements attribute
                """
                return self.elements

        def set_elements (self, elements) :
                """
                Sets the elements. Needed in write mode.

                @param elements : The element names of the atoms
                """
                self.elements = elements

        def set_name (self, name) :
                """
                Sets the name of the system. Needed in write mode
                """
                self.name = name

        def read_header (self) :
                """
                Set up info required for reading frames
                """
                line = self.file_object.readline()
                self.ntap = int(line.split()[0])
                if self.coords.shape == (0,3) :
                        self.coords = numpy.zeros((self.ntap,3)) 

                self.file_object.readline()

                elements = []
                for i in range(self.ntap) :
                        line = self.file_object.readline()
                        elements.append(line.split()[0])
                self.elements = elements

                self.file_object.seek(0)

        def get_plamsmol (self) :
                """
                Creates a PLAMS molecule object from the xyz-trajectory file
                """
                from scm.plams import Molecule
                oldposition = self.position
                self.rewind()
                coords, cell = self.read_next()
                plamsmol = Molecule.from_elements(self.elements)
                plamsmol.from_array(coords)

                # Return to original position
                self.rewind()
                for i in range(oldposition) :
                        self.read_next(read=False)

                return plamsmol

        def read_next (self, molecule=None, read=True) :
                """
                Reads the relevant info from current frame
                """
                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()

                cell = None
                # Read the coordinates
                self.read_coordinates(molecule)

                if self.firsttime :
                        self.firsttime = False

                self.position += 1
                
                return self.coords, cell

        def read_coordinates (self, molecule) :
                """
                Read the coordinates from file, and place them in the molecule
                """
                cell = None
                for i in range(2) :
                        line = self.file_object.readline()
                        if len(line) == 0 :
                                return None, None           # End of file is reached
                for i in range(self.ntap) :
                        line = self.file_object.readline()
                        self.coords[i,:] = [float(w) for w in line.split()[1:4]]

                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords, cell, molecule)

        def _is_endoffile (self) :
                """
                If the end of file is reached, return coords and cell as None
                """
                end = False
                for i in range(self.ntap+2) :
                        line = self.file_object.readline()
                        if len(line) == 0 :
                                end = True
                                break
                return end

        def write_next (self,coords=None,molecule=None,cell=[0.,0.,0.],energy=0.,step=None,conect=None) :
                """
                Write the coordinates to the next available spot in the file
                """
                if isinstance(molecule,Molecule) :
                        coords, cell, elements = self._read_plamsmol(molecule)[:3]
                        if self.position == 0 :
                                self.elements = elements
                cell = self._convert_cell(cell)

                self.write_moldata(coords, cell, energy, step)

                self.position += 1

        def write_moldata (self, coords, cell, energy, step) :
                """
                Write all molecular info to file
                """
                if step is None :
                        block = create_xyz_string(self.elements, coords)
                else :  
                        box = None
                        if cell is not None :
                                #box = PDBMolecule().box_from_vectors(cell)
                                box = cell_shape(cell)
                        block = create_xyz_string(self.elements,coords,energy,box,step,self.name)
                self.file_object.write(block)

        def _rewind_to_first_frame(self) :
                """
                Rewind the file to the first frame
                """
                self.file_object.seek(0)
                self.firsttime = True
                self.position = 0

        def _rewind_n_frames(self,nframes) :
                """
                Rewind the file by nframes frames
                """
                new_frame = self.position - nframes
                self._rewind_to_first_frame()
                for i in range(new_frame) :
                        self.read_next(read=False)


def create_xyz_string (elements, coords, energy=None, box=None, step=None, name='PlamsMol') :
        """
        Write an XYZ file based on the elements and the coordinates of the atoms
        """
        block = '%8i\n'%(len(elements))
        if step is not None :
                if energy is None : energy = 0.
                comment = '%-40s%6i %16.6f'%(name, step, energy)
                if box is not None :
                        for value in box :
                                comment += '%7.2f'%(value)
                block += comment
        block += '\n'
        for el,crd in zip(elements,coords) :
                block += '%8s '%(el)
                for x in crd :
                        block += '%20.10f '%(x)
                block += '\n'
        return block
