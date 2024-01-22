import numpy as np
import pandas as pd
a = 2.635
b = 1.755
d_HeC = 4.506
input = 'NPI2.pdb'
output = '3.pdb'
res_name = 'Mae'
C1,O1,O2 = 5,3,2
C2,O3,O4 = 6,1,4

def readpdb(pdb,res_name):
    inputfile = str(pdb)
    outputfile = inputfile.strip(".pdb")
    with open(inputfile,'r') as fp:
        content = fp.readlines()
        linesnumber = len(content)
    with open(outputfile+'.txt','w') as fp_w:
        lines = [] 
        for i in range (linesnumber):
            values = content[i].split()       
                # Extract values based on their positions in the format string
            if (values[0]=='ATOM' or values[0] == 'HETATM'):
                value1 = values[2] #atom_label
                value2 = res_name
                value3 = float(values[5]) #x
                value4 = float(values[6]) #y
                value5 = float(values[7]) #z
                value6 = values[10] #atom_note
                value7 = 1 #res_number
                newline = "%7s%7s%5d%8.3f%8.3f%8.3f%7s" % (
                    value1, value2, value7,value3, value4, value5, value6
                    )
                lines.append(newline+'\n')

        fp_w.writelines(lines)
    
    data = pd.read_csv(outputfile+'.txt',delim_whitespace=True,names=['Atom_label','Residue','Res_number','x','y','z','Note'])
    return data

def norm(v):
    l=np.linalg.norm(v)
    norm_v=v/l
    return norm_v

def pdbformat(label,xyz,res_name,count):
    x,y,z=xyz
    count = int(count)
    formatted_line = "%-6s%5d%5s%4s%10d%8.3f%8.3f%8.3f%6s%6s%7s" % (
                        'ATOM', count, label, res_name, 1, x, y, z,'1.0','0.0',label
            )  
    return formatted_line+'\n'

def pdbgetHe_Zr(a,b,d_HeC,pdb,num_C1,num_O1,num_O2,res_name,count):
    count = int(count)
    pd= readpdb(pdb,res_name)
    C1 = np.asarray(pd.loc[num_C1-1,['x','y','z']],dtype = float)
    O1 = np.asarray(pd.loc[num_O1-1,['x','y','z']],dtype = float)
    O2 = np.asarray(pd.loc[num_O2-1,['x','y','z']],dtype = float)

    X = 0.5*(O1+O2)
    C1X = X-C1
    norm_C1X = norm(C1X) #x axis, C1--He
    C1O1 = O1-C1
    C1Y = norm_C1X*np.dot(norm_C1X,C1O1.reshape(3,1))
    YO1=C1O1-C1Y  
    norm_YO1= norm(YO1) # y axis, vertical to C1--He

    He1 = C1+d_HeC*norm_C1X
    Zr1 = a*norm_C1X+b*norm_YO1+C1
    Zr2 = a*norm_C1X-b*norm_YO1+C1

    He_num = count+1
    Zr1_num = count+2
    Zr2_num = count+3

    He1_string = pdbformat("He",He1,res_name,He_num) 
    Zr1_string = pdbformat("Zr",Zr1,res_name,Zr1_num) 
    Zr2_string = pdbformat("Zr",Zr2,res_name,Zr2_num) 
    
    return(He1_string,Zr1_string,Zr2_string)

def get_linker_count(pdb):
    with open(pdb,'r') as fp:
        content = fp.readlines()
        lastline=content[-1].split() 
        count=int(lastline[1])
    return count

def pdbHe_Zr(a,b,d_HeC,input,output,C1,O1,O2,C2,O3,O4,res_name):
    count = get_linker_count(input)
    #print(count)
    He1_pdb,Zr1_pdb,Zr2_pdb = pdbgetHe_Zr(a,b,d_HeC,input,C1,O1,O2,res_name,count)
    newcount=count+3
    He2_pdb,Zr3_pdb,Zr4_pdb = pdbgetHe_Zr(a,b,d_HeC,input,C2,O3,O4,res_name,newcount)

    with open(input,'r') as f:
        original_lines = f.readlines()
        #original_lines = original_lines[1:-1] #remove last blank line

    with open(output,'w') as f:
        newpdb = []
        for i in original_lines:
            newpdb.append(i)
        f.writelines(newpdb)

    with open(output,'a') as f:
        f.write(He1_pdb)
        f.write(Zr1_pdb)
        f.write(Zr2_pdb)
        f.write(He2_pdb)
        f.write(Zr3_pdb)
        f.write(Zr4_pdb)


pdbHe_Zr(a,b,d_HeC,input,output,C1,O1,O2,C2,O3,O4,res_name)
