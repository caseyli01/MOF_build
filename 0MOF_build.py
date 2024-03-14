import MOF_build
import numpy as np
import pandas as pd
import datetime
import os
#import glob
startTime = datetime.datetime.now()

MLM_filename = MOF_build.build.readpdb("opt_linker.pdb")  
frame = MOF_build.build.readpdb("FM3M_primitive.pdb")
extra_termination = MOF_build.build.readpdb("methyl.pdb")
extra_point_index = 1
distance_extra_terimination = 1.7
#print(MLM_filename.shape)                                                  
defined_ATOM = 'He'
defined_ATOM_M = 'Zr' 
cutoff = 5  
cut_residue_name = 'CUT'
metal_residue_name = 'ZR1'
outfile_name = 'frame_MOF'
subfolder = 'frame_MOF'
path_residue = 'Residues'
 

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------#
newpath = os.path.abspath ('')+'/'+subfolder+'/'    # input file
residue_path = os.path.abspath ('')+'/'+path_residue+'/'
os.makedirs(newpath,exist_ok=True) 
os.makedirs(residue_path,exist_ok=True) 

points_n = MOF_build.build.findtop_frame(frame,defined_ATOM)
edge_length = []
for i in range(points_n):
    for j in range(points_n):
        if (i < j):
            edge = MOF_build.build.length_square(frame.loc[i],frame.loc[j])
            edge_length.append(edge)

frame_edge = min(edge_length)                                       
#frame_edge = 50   #should be the square of expected length

L_filename = MOF_build.input.getL_file(MLM_filename,defined_ATOM,defined_ATOM_M)
L_filename.to_csv(residue_path+'linker.txt', sep='\t', header = None)
M_filename = MOF_build.input.getM_file(MLM_filename,defined_ATOM)
Zr_linker = MOF_build.input.getZr_side_file(MLM_filename,defined_ATOM_M)
ATOM_index =MOF_build.build.findTOPinlinker(MLM_filename,defined_ATOM)
print('\n'+"linker atom numbers:   "+str(L_filename.shape[0]))
#print(M_filename,Zr_linker)


#find user-defined locating-assisted atom in linker
point_A = np.asarray(MLM_filename.loc[ATOM_index[0],['x','y','z']],dtype = float)
point_B = np.asarray(MLM_filename.loc[ATOM_index[1],['x','y','z']],dtype = float)

#get outer points positions 
#outer points xyz compared to point_A+compared to point_B --> new PMMP file                                                               
Cut_in_linker = MOF_build.build.find_points_in_cutoff(L_filename, point_A, point_B, cutoff).reset_index(drop=True)   
Extra = MOF_build.build.calculate_user_deifined_termination(M_filename,Cut_in_linker,extra_termination,extra_point_index,distance_extra_terimination)
PMMP_file = pd.concat([Cut_in_linker,Extra], ignore_index=True).reset_index(drop=True)
PMMP_file.to_csv(residue_path+'CUT.txt', sep='\t', header = None)
print('\n'+"(unit = A)Atoms in ligand in cutoff range around ligand top point: "+'\n',PMMP_file)

Metal_file = MOF_build.build.find_points_in_cutoff(MLM_filename, point_A, point_B, 4)
Metal_file = Metal_file[Metal_file['Atom_label'] == defined_ATOM_M].reset_index(drop=True)

#zoom  #frame atoms "He"
MM_l = MOF_build.build.length_square(M_filename.iloc[0,:],M_filename.iloc[1,:])
zoom = (MM_l/frame_edge)**0.5
df1 = MOF_build.build.get_zoomedframe(points_n,frame,zoom)
print('\n'+"edge length is:   "+str(round(frame_edge,2)))

linker_n = MOF_build.build.get_linker_number(edge_length,frame_edge)
print('\n'"number of linkers:   "+str(linker_n))

#residue_count = 1 remove He

#count =  1+points_n
count = 1
df_mof = MOF_build.build.calculate_MOF_linker(point_A,point_B,points_n,MM_l,df1,L_filename,Zr_linker,count)
df_mof.to_csv('example_MOF_in.gro', sep='\t', header = None, index = False)


# for Cut points residue_Term_count = ligand_n+1
Cut_count = linker_n+1
df_outpoints = MOF_build.build.calculate_MOF_CUT_TERM(point_A,point_B,points_n,MM_l,df1,PMMP_file,Metal_file,Zr_linker,Cut_count,cut_residue_name)
df_outpoints.to_csv('example_MOF_out.gro', sep='\t', header = None)
df_outpoints = df_outpoints.reset_index(drop = True)
#print(Cut_count)
if PMMP_file.shape[0] !=0:
    Metal_count = Cut_count+int(df_outpoints.shape[0]/PMMP_file.shape[0])
else:
    Metal_count = Cut_count
df_metal = MOF_build.build.get_framenode_octahedral(point_A,points_n,MM_l,df1,Metal_file,Metal_count,metal_residue_name,residue_path)

df_metal.to_csv('example_MOF_metal.gro', sep='\t', header = None)
df_metal = df_metal.reset_index(drop = True)

df_all = pd.concat([df_mof,df_outpoints,df_metal], ignore_index=True)
df_all.to_csv(newpath+outfile_name+'.txt', sep='\t', header = None, index = False)


MOF_build.output.outgro(df_all,newpath+outfile_name,points_n)
MOF_build.output.outxyz(newpath+outfile_name,points_n)
MOF_build.output.outpdb(newpath+outfile_name,points_n)

endTime = datetime.datetime.now()
print('\n'+"Time cost (s):   "+str(endTime-startTime))

MOF_build.output.clean()
MOF_build.residue.residues2xyz(residue_path)
MOF_build.residue.residues2pdb(residue_path)
MOF_build.residue.residues2gro(residue_path)
MOF_build.residue.clean(residue_path)
