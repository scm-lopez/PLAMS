#!/home4/tc/usr/bin/python

#############################################################################
#
# class pytrajectory
#
#############################################################################

import os
import struct
import array
import numpy
from scm.plams import Molecule, PlamsError
from .trajectoryfile import TrajectoryFile

__all__ = ['DCDTrajectoryFile']

class DCDTrajectoryFile (TrajectoryFile) :
        """
        class for representing a dcd trajectory file
        """

        def __init__ (self, filename=None, mode='rb', fileobject=None, ntap=None) :
                """
                Creates the DCDFile object

                @param ntap : Number of atoms (constant throughout trajectory)
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
                self.coords = numpy.zeros((self.ntap,3))

                # DCD specific attributes
                self.celldata = False                   # Whether the cell data is in the file (only set in read-mode)
                self.byte_order = '@'                   # native endian
                self.timesteps = (0,0,0)                # (nsteps,startstep,steps_between_saves)
                self.delta = 2.0                        # The timestep (if stored)
                self.namnf = 0                          # The number of free atoms (if stored)
                self.stepsize = 0                       # The number of bytes used to store a single conformation
                self.intsize = struct.calcsize('i')     # The number of bytes in an integer
                self.floatsize = struct.calcsize('f')   # The number of bytes in a float

                # Skip to the trajectory part of the file
                if self.mode == 'rb' :
                        self.read_header()
                elif self.mode == 'wb' :
                        self.write_header()

        def set_byteorder (self, byteorder) :
                """
                String. Can be either '@', '<', or '>', denoting native endian,
                        little endian or big endian.
                """
                self.byte_order = byteorder

        def read_variable (self,size) :
                """
                Reads a variable from a binary file

                @note: 
                        It assumes that the variable has been written on a machine
                        with the same bitsizes for the data
                """
                if size == 'i' :
                        varsize = self.intsize
                if size == 'f' :
                        varsize = self.floatsize
                else :
                        varsize = struct.calcsize(size)
                var_char = self.file_object.read(varsize)
                fmt = self.byte_order+size
                var = struct.unpack(fmt,var_char)[0]

                return var

        def write_variable (self,var,size) :

                var_char = struct.pack(size,var)
                self.file_object.write(var_char)

        def read_header (self) :

                # This should be an int32
                input_integer = self.read_variable('i')
                if input_integer != 84 :
                        print(input_integer)
                        print('BAD DCD FORMAT')
                        raise PlamsError

                # Read CORD string
                car4 = self.file_object.read(4).decode()
                if car4 != 'CORD' :
                        print('BAD DCD FORMAT')
                        raise PlamsError

                # Here we read the data of the timesteps
                # NSET: The number of sets of coordinates
                # ISTART: The starting timestep
                # NSAVC: The number of timesteps between dcd saves
                timesteps = array.array('i')
                timesteps.fromfile(self.file_object,3)
                self.timesteps = tuple(timesteps)
                #print 'timesteps: ',self.timesteps

                # 5 blanc integers
                dummies = array.array('i')
                dummies.fromfile(self.file_object,5)

                # NAMNF: The number of free atoms
                self.namnf = self.read_variable('i')

                # DELTA: The timestep
                delta = self.read_variable('f')

                # 10 blanc integers
                dummies = array.array('i')
                dummies.fromfile(self.file_object,10)

                # The end size of the first block
                input_integer = self.read_variable('i')
                if input_integer != 84 :
                        print('BAD DCD FORMAT')
                        raise PlamsError

                ############################################################

                # The size of the next block
                input_integer = self.read_variable('i')
        
                if (input_integer-4)%80 == 0 :
                        # NTITLE: The number of 80 char title strings
                        ntitle = self.read_variable('i')

                        for i in range(ntitle) :
                                car = self.file_object.read(80).decode()
                                #print car

                        # The end size of this block
                        input_integer = self.read_variable('i')
                else :
                        print('BAD DCD FORMAT')
                        raise PlamsError

                ############################################################

                # This should be an int32
                input_integer = self.read_variable('i')
                if input_integer != 4 :
                        print('BAD DCD FORMAT')
                        raise PlamsError

                # Read in the number of atoms
                self.ntap = self.read_variable('i')
                #print 'ntap: ',self.ntap
                if self.coords.shape == (0,3) :
                        self.coords = numpy.zeros((self.ntap,3))

                # This should be an int32
                input_integer = self.read_variable('i')
                if input_integer != 4 :
                        print('BAD DCD FORMAT')
                        raise PlamsError
                
                if self.namnf != 0 :
                        # This should be an int32
                        input_integer = self.read_variable('i')
                        if input_integer != (self.ntap-self.namnf)*4 :
                                print('BAD DCD FORMAT')
                                raise PlamsError
                
                        # This should be an int32
                        dummies = array.array('i')
                        dummies.fromfile(self.file_object,(self.ntap-self.namnf))
                        #print 'list',list(dummies)

                        # This should be an int32
                        input_integer = self.read_variable('i')
                        if input_integer != (self.ntap-self.namnf)*4 :
                                print('BAD DCD FORMAT')
                                raise PlamsError

                self.stepsize += 3 * self.ntap * self.floatsize
                self.stepsize += 5 * self.intsize

        def read_next (self,molecule=None,read=True) :
                """
                Reads the geometry at the cursor position
                """
                # Check for end of file
                # This should be an int32
                input_integer_char = self.file_object.read(self.intsize)
                if len(input_integer_char) > 0 :
                        if self.firsttime :
                                self.celldata = False
                                input_integer = struct.unpack(self.byte_order+'i',input_integer_char)[0]
                                if input_integer == 48 :
                                        self.celldata = True
                else :
                        return None, None

                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()
               
                # Read cell data 
                cell = numpy.zeros((3,3))
                if self.celldata :
                        cell_array = array.array('d')
                        try :
                                cell_array.fromfile(self.file_object,6)
                                if self.byte_order != '@' :
                                        cell_array.byteswap()
                        except EOFError :
                                return None, None
                        cell[numpy.tril_indices(3)] = cell_array
                        self.file_object.read(self.floatsize)

                        # Coords
                        if self.firsttime :
                                input_integer = self.read_variable('i')
                        else :
                                self.file_object.read(self.intsize)

                if self.firsttime :
                        if input_integer != 4*self.ntap :
                                print('BAD DCD FORMAT')
                                raise PlamsError
                
                xcoords = array.array('f')
                try :
                        xcoords.fromfile(self.file_object,self.ntap)
                        #xcoords = numpy.fromfile(self.file,dtype='f',count=self.ntap)
                        if self.byte_order != '@' :
                                xcoords.byteswap()
                except EOFError :
                        return None, None
              
                for i in range(2) : 
                        if self.firsttime :
                                input_integer = self.read_variable('i') 
                                if input_integer != 4*self.ntap :
                                        print('BAD DCD FORMAT')
                                        raise PlamsError
                        else :
                                self.file_object.read(self.intsize)

                ycoords = array.array('f')
                try:
                        ycoords.fromfile(self.file_object,self.ntap)
                        if self.byte_order != '@' :
                                ycoords.byteswap()
                except EOFError :
                        return None, None

                for i in range(2) :
                        if self.firsttime :
                                input_integer = self.read_variable('i')
                                if input_integer != 4*self.ntap :
                                        print('BAD DCD FORMAT')
                                        raise PlamsError
                        else :
                                self.file_object.read(self.intsize)

                zcoords = array.array('f')
                try :
                        zcoords.fromfile(self.file_object,self.ntap)
                        if self.byte_order != '@' :
                                zcoords.byteswap()
                except EOFError :
                        return None, None
              
                if self.firsttime : 
                        input_integer = self.read_variable('i') 
                        if input_integer != 4*self.ntap :
                                print('BAD DCD FORMAT')
                                raise PlamsError
                else :
                        self.file_object.read(self.intsize)

                # Check that self.coords has not been changed
                if not self.coords.shape == (self.ntap,3) :
                        raise PlamsError('The coordinate array has been changed outside the class')
                self.coords[:,0] = xcoords
                self.coords[:,1] = ycoords
                self.coords[:,2] = zcoords

                if self.firsttime :
                        if self.celldata :
                                self.stepsize += 6 * struct.calcsize('d')
                                self.stepsize += self.floatsize
                                self.stepsize += self.intsize
                        self.firsttime = False

                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords, cell, molecule)

                self.position += 1
              
                return self.coords, cell

        def _is_endoffile (self) :
                """
                If the end of file is reached, return coords and cell as None
                """
                test = self.file_object.read(self.stepsize)
                return len(test) < self.stepsize

        def write_header (self) :

                # This should be an int32
                self.write_variable(84,'i')

                # Read CORD string
                self.file_object.write('CORD'.encode())

                # Here we read the data of the timesteps
                # NSET: The number of sets of coordinates
                # ISTART: The starting timestep
                # NSAVC: The number of timesteps between dcd saves
                timesteps = array.array('i',[0,1,1])
                self.file_object.write(timesteps)

                # 5 blanc integers
                dummies = array.array('i',[0,0,0,0,0])
                self.file_object.write(dummies)

                # NAMNF: The number of free atoms
                self.write_variable(0,'i')

                # DELTA: The timestep
                delta = self.write_variable(self.delta,'f')

                # 10 blanc integers
                dummies = array.array('i')
                dummies.append(1)
                for i in range(8) :
                        dummies.append(0)
                dummies.append(24)
                self.file_object.write(dummies)

                # The end size of the first block
                self.write_variable(84,'i')

                ############################################################

                # The size of the next block
                self.write_variable(84,'i')
       
                ntitle = 1
                self.write_variable(ntitle,'i')
                title = 'REMARK Created with py_md'
                nlet = len(title)
                for i in range(80-nlet) :
                        title += ' '
                self.file_object.write(title.encode())

                self.write_variable(84,'i')
 
                ############################################################

                # This should be an int32
                self.write_variable(4,'i')

                # Read in the number of atoms
                self.write_variable(self.ntap,'i')

                # This should be an int32
                self.write_variable(4,'i')

        def write_next (self, coords=None, molecule=None, cell=[0.,0.,0.], energy=0.,step=None,conect=None) :
                """
                Writes the provided coordinates and cell vectors to the next DCD frame

                The coordinates and cell vectors can be provided directly, or in a molecule object
                """
                if isinstance(molecule,Molecule) :
                        coords, cell = self._read_plamsmol(molecule)[:2]
                cell = self._convert_cell(cell)

                if self.ntap != len(coords) :
                        print('BAD coordinate  list')
                        return

                # Invert the coords list [xcoords,ycoords,zcoords]
                #coords = numpy.array(coords).transpose()
                
                # This should be an int32
                self.write_variable(48,'i')
             
                cell = array.array('d',cell)
                self.file_object.write(cell)
                self.write_variable(0,'f')
                        
                # Coords
                self.write_variable(4*self.ntap,'i')
               
                xcoords = array.array('f',coords[:,0])
                self.file_object.write(xcoords)
                
                for i in range(2) : 
                        self.write_variable(4*self.ntap,'i')
               
                ycoords = array.array('f',coords[:,1])
                self.file_object.write(ycoords) 
                
                for i in range(2) :
                        self.write_variable(4*self.ntap,'i')
               
                zcoords = array.array('f',coords[:,2])
                self.file_object.write(zcoords) 

                self.write_variable(4*self.ntap,'i') 

                self.position += 1

        def _convert_cell (self, cell) :
                """
                Check format of cell, and convert to required lower triangle
                """
                if cell is None :
                        cell = numpy.zeros((3,3))[numpy.tril_indices(3)]
                elif len(cell)==3 and isinstance(cell[0],float) or isinstance(cell[0],numpy.float64) :
                        cell = [cell[0],0.,cell[1],0.,0.,cell[2]]
                else :
                        cell = numpy.array(cell)
                        cell = cell[numpy.tril_indices(3)]
                return cell

        def _rewind_to_first_frame(self) :
                """
                Rewind the file to the first frame
                """
                self.file_object.seek(0)
                self.firsttime = True
                self.position = 0
                self.stepsize = 0
                self.read_header()

        def _rewind_n_frames(self,nframes) :
                """
                Rewind the file by nframes frames
                """
                whence = os.SEEK_CUR
                stepsize = self.stepsize + self.intsize
                self.file_object.seek(-nframes*stepsize,whence)
                self.position = self.position - nframes
