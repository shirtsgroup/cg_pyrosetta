import cg_pyrosetta

def init(**kwargs):
    """
    Function used to initialize CG PyRosetta. This function calls the
    pyrosetta.init() function with the proper commandline values to 
    read in the CG models implemented in this package.
    """
    if kwargs is not None:
        for key in kwargs.keys():
            rosetta_key = "--" + key
            cg_pyrosetta.cmd_line_options_defaults[rosetta_key] = kwargs[key]
    
    init_string = ""
    for key in cg_pyrosetta.cmd_line_options_defaults.keys():
        if isinstance(cg_pyrosetta.cmd_line_options_defaults[key], list):
            for value in cg_pyrosetta.cmd_line_options_defaults[key]:
                init_string += key + " " + value + " "
        if isinstance(cg_pyrosetta.cmd_line_options_defaults[key], str):
            init_string += key + " " + cg_pyrosetta.cmd_line_options_defaults[key] + " "

    print(init_string)
    print(cg_pyrosetta.pyrosetta.__file__)
    cg_pyrosetta.pyrosetta.init(init_string)

