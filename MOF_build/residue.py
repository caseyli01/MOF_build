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
                    values = lines[i].split()
                    value_label = values[1] 
                    value_label = re.sub(r'\d', '', value_label)
                    try:                      
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

def clean(path):
    files = glob.glob(path+'*.txt')
    for f in files:
        os.remove(f)
    
