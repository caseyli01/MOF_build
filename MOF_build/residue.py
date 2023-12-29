import re
import glob
import os

def residues2xyz(path):
    files = glob.glob(path+'*.txt')
    for f in files:
        with open(f,'r') as fp:
            lines = fp.readlines()
            name = lines[0].split()[2]
            number = len(lines)
        with open(path+str(name)+'.xyz','w') as f:
            newxyz = []
            newxyz.append(str(number)+'\n'+name+'\n')
            for i in range (number):
                    # Split the line into individual values (assuming they are separated by spaces)
                    values = lines[i].split()

                    value_label = values[1] 
                    value_label = re.sub(r'\d', '', value_label)
                    try:
                        # Try to convert values[4] to a float
                        
                        value_x = float(values[4]) #x      
                        value_y = float(values[5]) #y
                        value_z = float(values[6]) #z
                        formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                                    value_label, value_x, value_y, value_z
                        )        
                        newxyz.append(formatted_line+'\n')        
                    except ValueError:
                        value_x = float(values[5]) #x      
                        value_y = float(values[6]) #y
                        value_z = float(values[7]) #z
                        formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                                    value_label, value_x, value_y, value_z
                        )        
                        newxyz.append(formatted_line+'\n')  

            f.writelines(newxyz)
            print(str(name)+'.xyz'+"    done")


def residues2pdb(path):
    files = glob.glob(path+'*.txt')
    for f in files:
        with open(f,'r') as fp:
            lines = fp.readlines()
            name = lines[0].split()[2]
            number = len(lines)
        with open(path+str(name)+'.pdb','w') as f:
            newpdb = []
            newpdb.append(str(number)+'\n'+name+'\n')
            for i in range (number):
                    # Split the line into individual values (assuming they are separated by spaces)
                    values = lines[i].split()

                    value_label = values[1] 
                    value_label = re.sub(r'\d', '', value_label)
                    try:
                        # Try to convert values[4] to a float
                        # Extract values based on their positions in the format string
                        value1 = 'ATOM'
                        value2 = i+1
                        value3 = values[1] #label
                        value4 = values[2] #residue
                        value5 = 1 #residue number
                        value6 = float(values[4]) #x      
                        value7 = float(values[5]) #y
                        value8 = float(values[6]) #z
                        value9 = '1.00'
                        value10 = '0.00'
                        value11 = values[6]
                        # Format the values using the specified format string
                        formatted_line = "%-6s%5d%5s%4s%10d%8.3f%8.3f%8.3f%6s%6s%7s" % (
                                    value1, value2, value3, value4, value5, value6, value7, value8,value9,value10,value11
                        )        
                        lines[i] = formatted_line+'\n'        
                        newpdb.append(formatted_line+'\n')        
                    except ValueError:
                        value1 = 'ATOM'
                        value2 = i+1
                        value3 = values[1] #label
                        value4 = values[2] #residue
                        value5 = 1 #residue number
                        value6 = float(values[5]) #x      
                        value7 = float(values[6]) #y
                        value8 = float(values[7]) #z
                        value9 = '1.00'
                        value10 = '0.00'
                        value11 = values[6]
                        # Format the values using the specified format string
                        formatted_line = "%-6s%5d%5s%4s%10d%8.3f%8.3f%8.3f%6s%6s%7s" % (
                                    value1, value2, value3, value4, value5, value6, value7, value8,value9,value10,value11
                        )        
                        lines[i] = formatted_line+'\n'        
                        newpdb.append(formatted_line+'\n')        


            f.writelines(newpdb)
            print(str(name)+'.pdb'+"    done")





def clean(path):
    files = glob.glob(path+'*.txt')
    for f in files:
        os.remove(f)
    
