import cg_pyrosetta

def init(**kwargs):
    if kwargs is not None:
        for key in kwargs.keys():
            cg_pyrosetta.cmd_line_options_defaults[key] = kwargs[key]
    
    init_string = ""
    for key in cg_pyrosetta.cmd_line_options_defaults.keys():
        if isinstance(cg_pyrosetta.cmd_line_options_defaults[key], list):
            for value in cg_pyrosetta.cmd_line_options_defaults[key]:
                init_string += key + " " + value + " "
        if isinstance(cg_pyrosetta.cmd_line_options_defaults[key], str):
            init_string += key + " " + cg_pyrosetta.cmd_line_options_defaults[key] + " "

    print(init_string)
    cg_pyrosetta.pyrosetta.init(init_string)

