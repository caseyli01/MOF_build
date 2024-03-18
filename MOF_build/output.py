import re
import glob
import os

def clean():
    files1 = glob.glob('example*')
    for f in files1:
        os.remove(f)

    files2 = glob.glob('*.txt')
    for f in files2:
        os.remove(f)


def get_box_dimension(file):
    x1,x2,y1,y2,z1,z2 = get_box(file)
    dx,dy,dz = abs(x1-x2), abs(y1-y2), abs(z1-z2)
    dimension = [str(round(dx,5)), str(round(dy,5)),str(round(dz,5))]
    s = ' '.join(dimension)
    return(s)

def get_box(file):
    x, y, z = [], [], []
    x, y, z = file['x'],file['y'],file['z']
    x1,y1,z1 = min(x), min(y), min(z)
    x2,y2,z2 = max(x), max(y), max(z)
    return(x1,x2,y1,y2,z1,z2)


def outgro(df_all,output,Hecount):
    with open(output+'.txt', 'r') as f:
        # Read the lines from the file
        lines = f.readlines()
        atoms_number = len(lines)

    newgro = []
    with open(output+'.gro', 'w') as fp:
        newgro.append("generated by MOF_BUILD"+'\n'+str(atoms_number-Hecount)+'\n')
        # Iterate over each line in the input file
        for i in range (Hecount,atoms_number):
            # Split the line into individual values (assuming they are separated by spaces)
            values = lines[i].split()
            # Extract values based on their positions in the format string
            #value1 = 'ATOM'
            value_atom_number = int(i+1-Hecount) #atom_number
            value_label = values[0] #atom_label
            value_resname = values[1] #residue_name
            value_resnumber = int(values[2]) #residue number
            value_x = float(values[3])/10 #x      
            value_y = float(values[4])/10 #y
            value_z = float(values[5])/10 #z
            #value11 = values[6] #note
            # Format the values using the specified format string
            formatted_line = "%5d%-5s%5s%5d%8.4f%8.4f%8.4f" % (
                        value_resnumber, value_resname, value_label, value_atom_number, value_x, value_y, value_z
            )        
            newgro.append(formatted_line+'\n')        

        tail = get_box_dimension(df_all.iloc[2:,3:6])+'\n'
        newgro.append(tail)
        fp.writelines(newgro)

def outxyz(output,Hecount):
    with open(output+'.txt', 'r') as f:
        # Read the lines from the file
        lines = f.readlines()
        atoms_number = len(lines)-Hecount

    newxyz = []
    with open(output+'.xyz', 'w') as fp:
        newxyz.append(str(atoms_number)+'\n'+"generated by MOF_BUILD"+'\n')
        # Iterate over each line in the input file
        for i in range (Hecount,atoms_number):
            # Split the line into individual values (assuming they are separated by spaces)
            values = lines[i].split()
            # Extract values based on their positions in the format string
            #value1 = 'ATOM'
            #value_atom_number = int(i+1) #atom_number
            value_label = values[0] #atom_label
            value_label = re.sub(r'\d', '', value_label)
            #value_resname = values[1] #residue_name
            #value_resnumber = int(values[2]) #residue number
            value_x = float(values[3]) #x      
            value_y = float(values[4]) #y
            value_z = float(values[5]) #z
            #value11 = values[6] #note
            # Format the values using the specified format string
            formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                        value_label, value_x, value_y, value_z
            )        
            newxyz.append(formatted_line+'\n')        

        fp.writelines(newxyz)

def outpdb(output,Hecount):
    with open(output+'.txt', 'r') as f:
    # Read the lines from the file
        lines = f.readlines()
        atoms_number = len(lines)

    newpdb = []
    with open(output+'.txt', 'w') as fp:
        # Iterate over each line in the input file
        for i in range (atoms_number):
            # Split the line into individual values (assuming they are separated by spaces)
            values = lines[i].split()
            # Extract values based on their positions in the format string
            value1 = 'ATOM'
            value2 = int(i+1-Hecount)
            value3 = values[0] #label
            value4 = values[1] #residue
            value5 = int(values[2]) #residue number
            value6 = float(values[3]) #x      
            value7 = float(values[4]) #y
            value8 = float(values[5]) #z
            value9 = '1.00'
            value10 = '0.00'
            value11 = values[6]
            # Format the values using the specified format string
            formatted_line = "%-6s%5d%5s%4s%10d%8.3f%8.3f%8.3f%6s%6s%7s" % (
                        value1, value2, value3, value4, value5, value6, value7, value8,value9,value10,value11
            )        
            lines[i] = formatted_line+'\n'        

        fp.writelines(lines)


    with open(output+'.pdb', 'w') as fp:
        # Iterate over each line in the input file
        newpdb.append("generated by MOF_BUILD"+'\n')
        newpdb.append(lines[Hecount])

        for i in range (Hecount+1,len(lines)):
            lastline = lines[i-1]
            thisline = lines[i]
            # Split the line into individual values (assuming they are separated by spaces)
            old_residue_number = lastline.split()[4]
            new_residue_number = thisline.split()[4]
            
            if(old_residue_number != new_residue_number):
                newline ='TER'+'\n'+thisline
                
            else:
                newline = thisline

            newpdb.append(newline)

        #print(newpdb_lines)
        fp.writelines(newpdb)
