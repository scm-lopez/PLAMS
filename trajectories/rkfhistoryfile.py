#!/usr/bin/env python

"""
Class representing an RKF MD file with changing numbers of atoms
"""

import numpy
from ..tools.periodic_table import PT
from ..mol.molecule import Molecule
from ..mol.atom import Atom
from ..core.settings import Settings
from .rkffile import RKFTrajectoryFile
from .rkffile import bohr_to_angstrom

__all__ = ['RKFHistoryFile']

class RKFHistoryFile (RKFTrajectoryFile) :

        def __init__ (self, filename, mode='rb', fileobject=None, ntap=None) :
                """
                Initializes the RKFHistoryFile object
                """
                self.added_atoms = None
                self.removed_atoms = None
                RKFTrajectoryFile.__init__(self,filename,mode,fileobject,ntap)

                self.input_elements = self.elements[:]
                self.versionhistory_length = 0
                self.system_versions = []
                self.version_history_items = []

        def set_elements (self, elements) :
                """
                Sets the elements. Needed in write mode.

                @param elements : The element names of the atoms
                """
                if self.position > 0 :
                        raise PLAMSError('Elements should not be changed while reading/writing is already in progress')
                RKFTrajectoryFile.set_elements(self, elements)
                self.input_elements = elements

        def read_header (self) :
                """
                Read the start molecule data from the InputMolecule section (not the Molecule section)
                """
                RKFTrajectoryFile.read_header (self, molecule_section='InputMolecule')

                # Now store the added and removed atoms along the trajectory
                # (This might be slow?)
                # I could also do it on the fly, but that may be messy when we move back and forth through the file
                if not 'SystemVersionHistory' in self.file_object.reader._sections :
                        return
                self.added_atoms = {}
                self.removed_atoms = {}
                version = 0
                for i in range(self.get_length()) :
                        new_version = self.file_object.read('History','SystemVersion(%i)'%(i+1))
                        if new_version == version :
                                continue
                        self.added_atoms[i] = {}
                        self.removed_atoms[i] = {}
                        # Now look for the added and removed atoms
                        removed_atoms = []
                        if 'RemovedAtoms(%i)'%(new_version) in self.file_object.reader._sections['SystemVersionHistory'] :
                                removed_atoms = self.file_object.read('SystemVersionHistory','RemovedAtoms(%i)'%(new_version))
                        added_atoms = []
                        if 'AddedAtoms(%i)'%(new_version) in self.file_object.reader._sections['SystemVersionHistory'] :
                                added_atoms = self.file_object.read('SystemVersionHistory','AddedAtoms(%i)'%(new_version))
                        # Now find the corresponding elements
                        chemSysNum = self.file_object.read('SystemVersionHistory','SectionNum(%i)'%(new_version))
                        sectionname = 'ChemicalSystem(%i)'%(chemSysNum-1)
                        if chemSysNum-1 == 0 :
                                sectionname = 'InputMolecule'
                        elements = [PT.get_symbol(atnum) for atnum in self.file_object.read(sectionname,'AtomicNumbers')]
                        elements = self.file_object.read(sectionname,'AtomSymbols')
                        removed_elements = [elements[i-1].split('.')[0] for i in removed_atoms]
                        sectionname = 'ChemicalSystem(%i)'%(chemSysNum)
                        elements = [PT.get_symbol(atnum) for atnum in self.file_object.read(sectionname,'AtomicNumbers')]
                        added_elements = [elements[i-1].split('.')[0] for i in added_atoms]
                        # Now store the elements
                        for iat,el in zip(removed_atoms,removed_elements) :
                                self.removed_atoms[i][iat] = el 
                        for iat,el in zip(added_atoms,added_elements) :
                                self.added_atoms[i][iat] = el
                        version = new_version

        # FIXME: The write_header section writes the starting molecule to the Molecule section, 
        #        not the final molecule (like the Fortran code does)
        def write_header (self, coords, cell) :
                """
                Write Molecule info to file (elements, periodicity)
                """
                # First write the general section
                self._write_general_section()
                
                # Then write the input molecule
                self._write_molecule_section(coords, cell)
                # I think the InputMolecule is mandatory in this case
                self._write_molecule_section(coords, cell, section='InputMolecule')
                if self.include_mddata :
                        # Start setting up the MDHistory section as well
                        self.file_object.write('MDHistory','blockSize',100)

        def read_coordinates (self, i, molecule, cell, bonds) :
                """
                Read the coordinates at step i
                """
                coords = self.coords.reshape(len(self.coords)*3)

                elements = self.read_elements_for_frame(i)
                if elements != self.elements :
                        self.elements = elements
                        self.coords = numpy.zeros((len(elements),3))
                        coords = self.coords.reshape((len(elements)*3))
                        # Rebuild the molecule (bonds will disappear for now)
                        if isinstance(molecule,Molecule) :
                                for at in reversed(molecule.atoms) :
                                        molecule.delete_atom(at)
                                molecule.properties = Settings()
                                for el in elements :
                                        atom = Atom(PT.get_atomic_number(el))
                                        molecule.add_atom(atom)

                coords[:] = self.file_object.read('History', 'Coords(%i)'%(i+1))
                coords *= bohr_to_angstrom
                # This changes self.coords behind the scenes

                # Assign the data to the molecule object
                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords,cell,molecule,bonds)

        def read_elements_for_frame (self, frame) :
                """
                Use the added and removed atoms to read the elements at each frame
                """
                elements = self.input_elements[:]
                for i in range(frame+1) :
                        if i in self.removed_atoms :
                                # Insert them into the elements list, in the correct place
                                elements = [el for iat,el in enumerate(elements) if not iat in self.removed_atoms[i]]
                        if i in self.added_atoms :
                                for iat,el in self.added_atoms[i].items() :
                                        elements.insert(iat,el)
                return elements

        def write_next (self, coords=None, molecule=None, elements=None, cell=[0.,0.,0.], conect=None, mddata=None) :
                """
                Write frame to next position in trajectory file

                mddata: A dictionary. It can contain the following keys:
                        'TotalEnergy', 'PotentialEnergy', 'Step', 'Velocities', 'KineticEnergy', 
                        'Charges', 'ConservedEnergy', 'Time', 'Temperature'
                """
                if isinstance(molecule,Molecule) :
                        coords, cell, elements, conect = self._read_plamsmol(molecule)
                # Make sure that the cell consists of three vectors
                cell = self._convert_cell(cell)
                if conect is not None :
                        if len(conect) == 0 : conect = None
                self.conect = conect

                # If this is the first step, write the header
                if self.position == 0 :
                        self.elements = elements
                        self.write_header(coords,cell)
                        # This is specific for the history-file object
                        self.input_elements = elements[:]

                # Define some local variables
                counter = 1
                step = self.position
                energy = 0.
                if mddata is not None :
                        if 'Step' in mddata :
                                step = mddata['Step']
                        if 'PotentialEnergy' in mddata :
                                energy = mddata['PotentialEnergy']

                # Write the history section
                counter = self.write_history_entry(step, coords, cell, conect, energy)

                if self.include_mddata and mddata is not None :
                        self.write_mdhistory_section(mddata)
                
                # If a change took place, write it.       
                if elements != self.elements or self.position==0 :
                        self.write_version_history(elements, coords, cell)
                        # Now update the elements
                        self.elements = elements
                        # Write the molecule sections
                        chemsysversion = len(self.system_versions)
                        self._write_molecule_section (coords, cell, section='ChemicalSystem(%i)'%(chemsysversion))
                        #self._write_molecule_section(coords, cell)
                counter = self.write_system_version_history_entry(counter)

                self.position += 1
                #if self.position%self.saving_freq == 0 : self.file_object.save()

        def write_version_history (self, elements, coords, cell) :
                """
                Write the version history
                """
                # Enter the correct SystemVersionHistory entry into the History section
                self.versionhistory_length += 1
                version = self.versionhistory_length

                # Find the corresponding ChemicalSystem and write it
                chemsysversion = self.check_for_chemical_system(elements)
                if chemsysversion is None :
                        self.system_versions.append(elements[:])
                        chemsysversion = len(self.system_versions)

                # Now add an entry to SystemVersionHistory
                added_atoms, removed_atoms = self.find_system_change(elements)
                self.file_object.write('SystemVersionHistory','nEntries',version)
                self.file_object.write('SystemVersionHistory','currentEntryOpen',False)
                if 'SectionNum' not in self.version_history_items :
                        self.version_history_items.append('SectionNum')
                svh_counter = self.version_history_items.index('SectionNum') + 1
                self._write_keydata_in_history('SectionNum',svh_counter,False,1,version,chemsysversion,'SystemVersionHistory')
                # Avoid writing empty data into this section by storing the items.
                if len(added_atoms) > 0 :
                        if 'AddedAtoms' not in self.version_history_items :
                                self.version_history_items.append('AddedAtoms')
                        svh_counter = self.version_history_items.index('AddedAtoms') + 1
                        data = [iat+1 for iat in added_atoms]
                        self._write_keydata_in_history('AddedAtoms',svh_counter,False,len(data),version,data,'SystemVersionHistory')
                if len(removed_atoms) > 0 :
                        if 'RemovedAtoms' not in self.version_history_items :
                                self.version_history_items.append('RemovedAtoms')
                        svh_counter = self.version_history_items.index('RemovedAtoms') + 1
                        data = [iat+1 for iat in removed_atoms]
                        self._write_keydata_in_history('RemovedAtoms',svh_counter,False,len(data),version,data,'SystemVersionHistory')

        def write_system_version_history_entry (self, counter) :
                """
                Write the entry for the SystemVersionHistory into the History section
                """
                version = self.versionhistory_length
                self._write_keydata_in_history('SystemVersion', counter, False, 1, self.position+1, version)
                counter += 1
                return counter

        def find_system_change (self, elements) :
                """
                Find out which atoms were added and/or deleted
                """
                # First find out which elements were removed
                removed_atoms = []
                position = 0
                for i,el in enumerate(self.elements) :
                        if el == elements[position] :
                                position += 1
                        else :
                                removed_atoms.append(i)

                # Then find out which elements were added
                added_atoms = [i for i in range(position,len(elements))]
                
                # Now store them (this is actually not really necessary)
                #for iat in removed_atoms :
                #        self.removed_atoms[self.position][iat] = self.elements[iat]
                #for iat in added_atoms :
                #        self.added_atoms[self.position][iat] = elements[iat]
                return added_atoms, removed_atoms

        def check_for_chemical_system (self, elements) :
                """
                Check if the new chemical system was encountered before
                """
                version = None
                for i,prev_elements in enumerate(self.system_versions) :
                        if elements == prev_elements :
                                version = i+1
                return version

