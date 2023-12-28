# MOF_build

MOF_build package with MOF_build module, can generate MOF with linker, suitable to build defective model, with user-defined cutoff termination, and metal cluster.(UiO-66 or other simlilar octahedral cluster center)
This version corrects the Zr cluster position based on the framework structure.

Lines in script need to be modified are as below:


--- 

MLM_filename = MOF_build.build.readpdb("opt_linker.pdb") \
frame = MOF_build.build.readpdb("FM3M_primitive3.pdb") \
extra_termination = MOF_build.build.readpdb("methyl.pdb") \
extra_point_index = 1 \
distance_extra_terimination = 1.7 \
defined_ATOM = 'He' \
defined_ATOM_M = 'Zr' \
cutoff = 5 \
cut_residue_name = 'CUT' \
metal_residue_name = 'ZR1' \
outfile_name = 'frame_MOF' \
subfolder = 'frame_MOF' \
path_residue = 'Residues' 

--- 


Need a linker file (MLM_filename) with a pair of location-assisted atoms, the location-assisted atoms form the framework structure (frame),
which only need the relative struture, whose edge length don't need to be equal to the real edge length in MOF.

>MLM_filename = MOF_build.build.readpdb("linker.pdb") \
  frame = MOF_build.build.readpdb("FM3M_primitive3.pdb")      


**The location-assisted atoms in the prepared *.pdb file have to be with the same atom name (atom label), and have to be unique enough in case to make the programme gets lost in the huge candidates.**\
**NEVER USE common atoms, like C O N H etc.**

The connector metal cluster atom_name of MOF should be announced in the start to separate the linker and metal-cluster for further generation.

>defined_ATOM = 'He'\
  defined_ATOM_M = 'Zr'

Set a cutoff value for the un-saturated points, like top points. it will only cut the linker file without the metal(like 'Zr1')

>cutoff = 5


Define your preferred name for the cut termination residues and metal cluster residues.Three letter name is adviced in case some bugs with the output file format havent been fixed(especially output.pdb)

>cut_residue_name = 'CUT'\
  metal_residue_name = 'MET'

Add you preferred termination to the cutoff and set the distance to put the target atom(and the whole termination) away from your cutoff.

>extra_termination = MOF_build.build.readpdb("methyl.pdb")\
>extra_point_index = 1 #first atom in the termination *.pdb file\
>distance_extra_terimination = 1.7

output files like *.pdb, *.gro and *.xyz will be stored in the subfolder and named by the outfile_name.

>outfile_name = 'frame_MOF'\
>subfolder = 'frame_MOF'

Residues appeared in the final structure will be stored in a subfolder whose name is "path_residue = 'Residues'" for further forcefield generation.\
**better to optimize the residue addtionally, especially for the cutoff.**

>path_residue = 'Residues'
