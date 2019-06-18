import shutil
import os


class PyRosettaBuilder():
    """
    PyRosettaBuilder will be an object used build and edit a CG version of
    pyrosetta.
    """
    
    def __init__(self, clean_pyrosetta_path, pyrosetta_path, inputs):
        if not os.path.isdir(pyrosetta_path):
            self.pyrosetta_path = self.copyPyRosetta(clean_pyrosetta_path, pyrosetta_path)
        else:
            self.pyrosetta_path = os.path.abspath(pyrosetta_path)
        self.inputs = inputs


    def copyPyRosetta(self, original, copy):  # VERY SLOW
        """
        Copy's a working pyrosetta directory into a new file system
        such that we don't ruin a working copy of PyRosetta
        """
        shutil.copytree(original, copy)
        return(os.path.abspath(copy))


    def buildCGPyRosetta(self):
        """
        function which takes in the path to a vanila PyRosetta package and
        adds CG functionality via copying files
        
        Arguments
        ---------

        path : path to the pyrosetta package being edited
        inputs : input directory with specific file structure denoting residue types and atom types
        """

        # Creating files where new residues will be injected

        if not os.path.exists(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','cg_models')):
            os.mkdir(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','cg_models'))

        if not os.path.exists(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','cg_models')):
            os.mkdir(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','cg_models'))

        if not os.path.exists(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models')):
            os.mkdir(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models'))

        self.addAtomTypes(os.path.join(self.inputs, 'atom_type_sets'), header=True)
        self.turnOffExtras()
        self.addResidueTypesNew(os.path.join(self.inputs, 'residue_type_sets'), header=True)
        self.addMMTorsionTypes(os.path.join(self.inputs, 'mm_atom_type_sets'))
        self.addPatches(os.path.join(self.inputs, 'residue_type_sets', 'patches'), header = True)
        self.addMMAtomTypes(os.path.join(self.inputs, 'mm_atom_type_sets'))
        self.addMMAngleTypes(os.path.join(self.inputs, 'mm_atom_type_sets'))



    def unBuildCGPyRosetta(self, path):
        """
        Intentions is to 'unbuild' our CG PyRosetta
        """
        pass



    def addAtomTypes(self, path, header = False):
        """
        Method for adding atoms from a specific path to the modified PyRosetta filesystem, checks for duplicates

        Arguments
        ---------

        self: PyRosettaBuilder Class
        path: str path of where input atom_properties.txt
        header: Bool used to determine whether or not write the header showing where custom atoms are added
        """
        input_path = os.path.abspath(path)
        
        atom_lines = []
        for files in os.listdir(path):
            if files.endswith(".txt"):
                with open(os.path.join(input_path,files), 'r') as f:  
                    for line in f.readlines()[1:]: # Skip header of atom_properties.txt files
                        atom_lines.append(line)


        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','atom_properties.txt'), 'r') as atom_file:
            original_file = atom_file.readlines()
            original_file = original_file[:172]
            
        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','atom_properties.txt'), 'w') as atom_file: 
            atom_file.writelines(original_file)

        # opening atom_properties.txt and appending new atom lines
        if header:
            with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','atom_properties.txt'), 'a') as atom_file:
                atom_file.write('\n')
                atom_file.write('##############################\n')
                atom_file.write('## Custom Added Atom Types ###\n')
                atom_file.write('##############################\n')
                atom_file.write('\n')
                    

        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','atom_properties.txt'), 'a') as atom_file:
            for atom_line in atom_lines:
                atom_file.write(atom_line)

    
    def addMMTorsionTypes(self, path):
        """
        Method for adding torsion definitions from a specific path to the modified PyRosetta filesystem, checks for duplicates

        Arguments
        ---------

        self: PyRosettaBuilder Class
        path: str path of where input mm_torsion_params.txt
        header: Bool used to determine whether or not write the header showing where custom atoms are added
        """
        input_path = os.path.abspath(path)
        
        torsion_lines = []

        with open(os.path.join(input_path,'mm_torsion_params.txt'), 'r') as f:  
            for line in f.readlines()[1:]: # Skip header of atom_properties.txt files
                torsion_lines.append(line)
        
        # Write lines to modified pyrosetta/.../mm_torsion_params.txt
        # with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','mm_torsion_params.txt'), 'r') as torsion_file:
        #    previous_lines = torsion_file.readlines()


        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','mm_torsion_params.txt'), 'w') as torsion_file:
            for torsion_line in torsion_lines:
                # Ensure there are no duplicate atom_properties lines
                    torsion_file.write(torsion_line)

    def addMMAtomTypes(self, path):
        """
        Method for adding MM AtomTypes into PyRosetta, required for evaluating mm_lj_atr/rep/etc. and for use in mm_twist and mm_bend terms

        Arguments
        ---------

        self: PyRosettaBuilder Class
        path: str path of where input mm_torsion_params.txt
        header: Bool used to determine whether or not write the header showing where custom atoms are added
        """
        input_path = os.path.abspath(path)
        
        mm_atom_lines = []

        # Reading each line of 'mm_atom_properties.txt'
        with open(os.path.join(input_path,'mm_atom_properties.txt'), 'r') as f:  
            for line in f.readlines()[1:]: # Skip header of atom_properties.txt files
                mm_atom_lines.append(line)
        
        # Write lines to modified pyrosetta/.../mm_torsion_params.txt
        #with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','mm_atom_properties.txt'), 'r') as mm_atom_file:
        #    previous_lines = mm_atom_file.readlines()


        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','mm_atom_properties.txt'), 'w') as mm_atom_file:
            for mm_atom_line in mm_atom_lines:
                mm_atom_file.write(mm_atom_line)

    def addMMAngleTypes(self, path):
            """
            Method for adding bond angle definitions from a specific path to the modified PyRosetta filesystem, checks for duplicates

            Arguments
            ---------

            self: PyRosettaBuilder Class
            path: str path of where input mm_torsion_params.txt
            header: Bool used to determine whether or not write the header showing where custom atoms are added
            """
            input_path = os.path.abspath(path)
            
            angle_lines = []
            
            with open(os.path.join(input_path,'mm_angle_params.txt'), 'r') as f:  
                for line in f.readlines(): # Skip header of atom_properties.txt files
                    angle_lines.append(line)
            
            
            # Write lines to modified pyrosetta/.../mm_torsion_params.txt
            # with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','par_all27_prot_na.prm'), 'r') as angle_file:
            #    previous_lines = angle_file.readlines()
            
            # start_angles = previous_lines.index('ANGLES\n')
            # print(start_angles)
            # for angle_line in angle_lines:
            #     if angle_line not in previous_lines:
            #         previous_lines.insert(start_angles+1, angle_line)
            #     else:
            #         print('Skipping MM Angle:', angle_line)

            with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','mm_atom_type_sets','cg_models','par_all27_prot_na.prm'), 'w') as angle_file:
                angle_file.writelines(angle_lines)
     
    def turnOffExtras(self):
        """
        Method used for commenting the extras.txt file in pyrosetta.

        Checks to see if file is already commented, if so gives a warning
        """
        # comment out extras.txt files

        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','extras.txt'), 'r') as exf:
            extras = exf.readlines()
        if extras[0][0]  == '#':
            print("PyRosetta's extra.txt file is already off!")
        else:
            extras = ['# '+line for line in extras]
            with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','extras.txt'), 'w') as exf:
                exf.writelines(extras)
    
    def turnOnExtras(self):
        """
        Method used for un-commenting the extras.txt file in pyrosetta.

        Checks to see if file is already un-commented, if so gives a warning
        """
        # comment out extras.txt files

        with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','extras.txt'), 'r') as exf:
            extras = exf.readlines()
        if extras[0][0]  != '#':
            print("PyRosetta's extra.txt file is already on!")
        else:
            extras = [line[2:] for line in extras]
            with open(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','extras.txt'), 'w') as exf:
                exf.writelines(extras)

    def addResidueTypes(self, path, header = False):
        """
        Method used for adding residue.param files to a modified PyRosetta4 filesystem

        Arguments
        ---------

        path : str path of where residue.param files are located
        header : Bool to write a header for the custom residue types locations
        """
        # add residue types into cg_models residue types
        custom_residue_path = os.path.abspath(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','cg_models','residue_types','custom'))

        if not os.path.isdir(custom_residue_path):
            os.mkdirs(os.path.join(self.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','cg_models','residue_types','custom'))
        
        for files in os.listdir(path):
            if files.endswith('.params'):
                if files in os.listdir(custom_residue_path):
                    print('Updating Residue:', files)
                shutil.copy(os.path.join(path, files), custom_residue_path)
                

        # add residue types to 'residue_types.txt' at the above l-caa residues (give priority for io strings)
        
        with open(os.path.join(self.pyrosetta_path, 'pyrosetta','database','chemical','residue_type_sets','cg_models','residue_types.txt'), 'r') as rtf:
            residue_type_lines = rtf.readlines()

        if header:
            residue_type_lines.insert(6,'### custom residues\n')
        # writting custom lines first
        for files in os.listdir(path):
            if os.path.join('residue_types','custom',files)+'\n' not in residue_type_lines and files.endswith('.params'):
                residue_type_lines.insert(7, os.path.join('residue_types','custom',files)+'\n')
            else:
                print('Skipping Residue: ', os.path.join('residue_types','custom',files))

        # Rewrite the residue_types.txt file with modifications
        with open(os.path.join(self.pyrosetta_path, 'pyrosetta','database','chemical','residue_type_sets','cg_models','residue_types.txt'), 'w') as rtf:
            rtf.writelines(residue_type_lines)
    
    def addResidueTypesNew(self, path, header = False):
        """
        Method for adding new residue types to a PyRosetta4 build

        Parameters
        ----------
        self: PyRosettaBuilder Class
            object running addResidueTypes
        path : str
            string indicating path to direcctory where residue.param files are
        header : Bool
            whether to include a header denoting custom residues
 
        """

        rel_path = os.path.relpath(path, os.path.join(self.pyrosetta_path, 'pyrosetta', 'database', 'chemical', 'residue_type_sets','cg_models'))

        with open(os.path.join(self.pyrosetta_path, 'pyrosetta','database','chemical','residue_type_sets','cg_models','residue_types.txt'), 'w') as rtf:

            if header:
                rtf.write('### custom residues\n')
            # writting custom lines first
            for residue in os.listdir(path):
                rtf.write(os.path.join(rel_path, residue)+'\n')


    def addPatches(self, path, header = False):
        """
        Method for adding new patches to a PyRosetta4 build

        Parameters
        ----------
        self : PyRosettaBuilder Class
            object running addPatches
        path : str
            string indicating dir where patch.txt files are
        header : Bool
            whether to include a header denoting custom patches

        """
        rel_path = os.path.relpath(path, os.path.join(self.pyrosetta_path, 'pyrosetta', 'database', 'chemical', 'residue_type_sets','cg_models'))

        with open(os.path.join(self.pyrosetta_path, 'pyrosetta','database','chemical','residue_type_sets','cg_models','patches.txt'), 'w') as rtf:
        
            if header:
                rtf.write('### custom patches residues\n')
            # writting custom lines first
            for patch in os.listdir(path):
                rtf.write(os.path.join(rel_path, patch)+'\n')
            


def add_import_path(path):
    """
    singler liner to import custom version of pyrosetta given any path
    """
    import sys
    sys.path.insert(0, path)


def main():
    pyrosetta_path = os.path.abspath('PyRosetta4.modified')
    input_path = 'inputs'
    
    builder = PyRosettaBuilder(pyrosetta_path, input_path)
    builder.buildCGPyRosetta()

if __name__ == '__main__':
    main()