def getL_file(MLM_filename,defined_ATOM,defined_ATOM_M):
    L_filename = MLM_filename[MLM_filename['Atom_label'] != defined_ATOM]
    L_filename = L_filename[L_filename['Atom_label'] != defined_ATOM_M]
    L_filename = L_filename.reset_index(drop=True)
    return L_filename

def getM_file(MLM_filename,defined_ATOM):
    M_filename = MLM_filename[MLM_filename['Atom_label'] == defined_ATOM]
    return M_filename

def getZr_side_file(MLM_filename,defined_ATOM_M):
    Zr_file = MLM_filename[MLM_filename['Atom_label'] == defined_ATOM_M].reset_index(drop=True)
    Zr_onepair = Zr_file.head(2)
    return Zr_onepair