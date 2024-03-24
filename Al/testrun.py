import numpy as np
import pandas as pd
import quaternion
import re
import datetime
startTime = datetime.datetime.now()



np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)

def get_box(file):
    x, y, z = [], [], []
    x, y, z = file['x'],file['y'],file['z']
    x1,y1,z1 = min(x), min(y), min(z)
    x2,y2,z2 = max(x), max(y), max(z)
    return(x1,x2,y1,y2,z1,z2)

def length_square(A1,A2):
    A1,A2 = list(A1),list(A2)
    x1,x2,y1,y2,z1,z2 = float(A1[3]),float(A2[3]),float(A1[4]),float(A2[4]),float(A1[5]),float(A2[5])
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    length2 = dx**2+dy**2+dz**2
    return length2

def readpdb(pdb):
    inputfile = str(pdb)
    outputfile = inputfile.strip(".pdb")
    with open(inputfile,'r') as fp:
        content = fp.readlines()
        linesnumber = len(content)
    lines = [] 
    with open(outputfile+'.txt','w') as fp_w:
        for i in range (linesnumber):
            # Split the line into individual values (assuming they are separated by spaces)
            values = content[i].split()       
            # Extract values based on their positions in the format string
            if (values[0]=='ATOM' or values[0] == 'HETATM'):
                value1 = values[2] #atom_label
                value2 = values[3] #res_name
                value3 = float(values[5]) #x
                value4 = float(values[6]) #y
                value5 = float(values[7]) #z
                value6 = values[10] #atom_note
                value7 = int(values[4])
            # Format the values using the specified format string
                newline = "%7s%7s%5d%8.3f%8.3f%8.3f%7s" % (
                    value1, value2, value7,value3, value4, value5, value6
                    )
     
                lines.append(newline+'\n')
        fp_w.writelines(lines)
    data = pd.read_csv(outputfile+'.txt',sep='\s+',names=['Atom_label','Residue','Res_number','x','y','z','Note'])
    return data

def normalize_vector(v):
    norm_v=v/np.linalg.norm(v)
    return norm_v

def findTOPinlinker(df,defined_atom_type):
    matching_indices = df[df['Atom_label'] == defined_atom_type].index.tolist()
    count_ATOMS = len(matching_indices)
    print('\n'+str(count_ATOMS)+' top points in ligand file')
    return matching_indices

def findtop_frame(df,defined_ATOM):
    matching_indices = df[df['Atom_label'] == defined_ATOM].index.tolist()
    count_ATOMS = len(matching_indices)
    print('\n'+str(count_ATOMS)+' top points in framework')
    return count_ATOMS

def get_linker_number(list,length):
    count = 0
    for i in list:
        if (length-1) < i <(length+1):
            count += 1
    return count

def find_points_in_radius(points, center, radius):
    points_identity = points.iloc[:, [0,1,2,5]]
    points_positions = np.array(points.loc[:, ['x','y','z']])
    center = np.array(center)
    distances = np.linalg.norm(points_positions - center, axis=1)
    indices = np.where(distances <= radius)[0]
    return points[indices]

def search_unique_vector(df1,MM_l):
    vector = []
    points_n = df1.shape[0]
    for i in range(points_n):            
        for j in range(points_n):
            x = round(length_square(df1.loc[i],df1.loc[j]))
            if ( x == round(MM_l)):
                point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
                point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
                v = point_B_frame -point_A_frame
                vAB = np.round(v, 1)
                vector.append(vAB)
    unique_vector = np.unique(vector, axis=0, return_index=False)
    unique_vector = pd.DataFrame(unique_vector,columns=['vx','vy','vz'])
    return unique_vector

def search_neighbor_vector(df1,i,MM_l):
    neighbor_vector = []
    points_n = df1.shape[0]         
    for j in range(points_n):
        x = round(length_square(df1.loc[i],df1.loc[j]))
        if ( x == round(MM_l)):
            point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
            point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
            v = point_B_frame -point_A_frame
            vAB = np.round(v, 1)
            neighbor_vector.append(vAB)
    neighbor_vector = pd.DataFrame(neighbor_vector,columns=['vx','vy','vz'])
    return neighbor_vector

def filtered_term_vector(df1,MM_l):
    points_n = df1.shape[0] 
    term_vector = []
    true_indices = []
    for i in range(points_n):
        unique_vector = search_unique_vector(df1,MM_l)
        neighbor_vector = search_neighbor_vector(df1,i,MM_l)
        df_vector = pd.concat([unique_vector,neighbor_vector],ignore_index=False)
        filtered_vector = df_vector.drop_duplicates(keep=False)
        term_vector.append(filtered_vector.to_numpy())
        true_indices.append(filtered_vector.index)
    return [term_vector,true_indices]

def calculate_q_rotation_with_axis_degree(axis,theta): #axis is HE---HE ,theta from O1--AXIS--O1'
    w = theta/2
    s = np.sin(w)
    q_real= np.array([np.cos(w)])
    q_ijk = s*axis
    q_r = np.concatenate([q_real,q_ijk])
    q_r = quaternion.from_float_array(q_r)
    return q_r

def calculate_q_rotation_with_vectors(p1,p2):
    q1 = quaternion.from_vector_part(p1)
    q2 = quaternion.from_vector_part(p2)
    r = q2*q1.conjugate()
    return r

def calculate_angle_rad(axis,p1, p2):
    axis = normalize_vector(axis)
    a_square=np.linalg.norm(p1)*np.linalg.norm(p1)-np.dot(p1,axis)*np.dot(p1,axis)
    b_square=np.linalg.norm(p2)*np.linalg.norm(p2)-np.dot(p2,axis)*np.dot(p2,axis)
    c_square=np.linalg.norm(p2-p1)*np.linalg.norm(p2-p1)
    a,b = np.sqrt(a_square),np.sqrt(b_square)
    cos_theta = (a_square+b_square-c_square)/(2*a*b)
    print(cos_theta,a_square,b_square,c_square,axis,p1,p2)
    cos_theta = np.clip(cos_theta,-1,1)
    theta_rad = np.arccos(cos_theta)
    print(theta_rad)
    return theta_rad

def points_generator(x_num,y_num,z_num,dx_value,dy_value,dz_value): 
    '''this function is to generate a group of 3d SCATTER defined by user for further grouping points'''
    dx = dx_value*np.array([[1,0,0]]) #dx_value works as a scalar
    dy = dy_value*np.array([[0,1,0]])
    dz = dz_value*np.array([[0,0,1]])
    # add x layer
    points = np.array([[0,0,0]])
    for i in range(1,x_num):
        points = np.concatenate((points,i*dx),axis=0)
    # add y layer
    points_x =points
    for i in range(1,y_num):
        points = np.concatenate((points,points_x+i*dy),axis = 0)
    # add z layer 
    points_xy = points
    for i in range(1,z_num):
        points = np.concatenate((points,points_xy+i*dz),axis = 0)
    
    return points

def groupA_one_step_to_groupB(first_B,group_A,d,points):
    point_dx = group_A+d
    group_B = first_B
    for i in range(group_A.shape[0]):
        point = np.reshape(point_dx[i], (1, 3))
        group_B = np.concatenate((group_B,point),axis = 0) if np.all(points == point, axis=1).any() else group_B
        group_B = np.unique(group_B,axis=0)
    return group_B

def groupA_one_step_to_groupA(group_A,d,points):
    point_dx = group_A+d
    for i in range(group_A.shape[0]):
        point = np.reshape(point_dx[i], (1, 3))
        group_A = np.concatenate((group_A,point),axis = 0) if np.all(points == point, axis=1).any() else group_A
        group_A = np.unique(group_A,axis=0)
    return group_A

def group_points_AB(x_num,y_num,z_num,dx_value,dy_value,dz_value):
   points=points_generator(x_num,y_num,z_num,dx_value,dy_value,dz_value)
   O = np.array([[0,0,0]])
   dx = dx_value*np.array([[1,0,0]]) #dx_value works as a scalar
   dy = dy_value*np.array([[0,1,0]])
   dz = dz_value*np.array([[0,0,1]])

   first_A = O
   first_B = O+dy
   group_A = first_A
   group_B = first_B

   for i in range(points.shape[0]):
      if group_A.shape[0]+group_B.shape[0] == points.shape[0]:
         break
      else:
         group_B_dx = groupA_one_step_to_groupB(first_B,group_A,dx,points)
         group_B_dy = groupA_one_step_to_groupB(first_B,group_A,dy,points)
         group_B_dz = groupA_one_step_to_groupA(group_B,dz,points)
         group_B_dxyz = [group_B,group_B_dx,group_B_dy,group_B_dz]
         group_B = np.concatenate([arr for arr in group_B_dxyz if arr is not None],axis = 0) 
         group_B = np.unique(group_B,axis=0)
         group_A_dx = groupA_one_step_to_groupB(first_A,group_B,dx,points)
         group_A_dy = groupA_one_step_to_groupB(first_A,group_B,dy,points)
         group_A_dz = groupA_one_step_to_groupA(group_A,dz,points)
         group_A_dxyz = [group_A,group_A_dx,group_A_dy,group_A_dz]
         group_A = np.concatenate([arr for arr in group_A_dxyz if arr is not None],axis = 0)
         group_A = np.unique(group_A,axis=0)
   
   return group_A,group_B

def get_center_point_of_face(p1_face,p2_face,p3_face):
    center_point = (normalize_vector(p1_face)+
                              normalize_vector(p2_face)+
                              normalize_vector(p2_face))/3
    return center_point

def find_solution(pAl1,pAl2,pAl1_1,pAl1_2,pAl1_3):
    Al1_Al2 = pAl2-pAl1
    vAl1_Al2 = normalize_vector(Al1_Al2)
    v12_1,v12_2,v12_3  = pAl1_1-pAl1,pAl1_2-pAl1,pAl1_3-pAl1
    v12_1,v12_2,v12_3= normalize_vector(v12_1),normalize_vector(v12_2),normalize_vector(v12_3)                        
    arr_1_2=np.vstack((v12_1,v12_2,v12_3))
    arr_1_2 = arr_1_2.astype(np.float64)
    vAl1_Al2 = vAl1_Al2.astype(np.float64)
    solution_1_2=np.dot(vAl1_Al2,np.linalg.inv(arr_1_2))
    return solution_1_2,arr_1_2

def get_rotate_array(arr,q):
    q_arr= quaternion.from_vector_part(arr)
    rotated_q_arr = q*q_arr*q.inverse()
    rotated_arr = quaternion.as_vector_part(rotated_q_arr)
    return rotated_arr

def get_axis2(solution_1_2,arr_1_2,solution_1_3,arr_1_3):
    axis1 = np.dot(solution_1_2,arr_1_2)
    axis2 = np.dot(solution_1_3,arr_1_3)
    q_axis = calculate_q_rotation_with_vectors(axis1,axis2)
    dx = np.array([1,0,0])
    axis0=quaternion.from_vector_part(dx)
    new_axis = q_axis*axis0
    new_axis_vector = quaternion.as_vector_part(new_axis)
    print(new_axis_vector)
    return new_axis_vector
 

def calculate_node(Metal_file,linker_cut_count,Residue_name,group_A,group_B,new_node_A,new_node_B):
#rotate as group, translate as group 
    Metal_count = linker_cut_count
    df_node = pd.DataFrame()
    for i in group_A:
        new_positions=new_node_A+i
        df_left = pd.DataFrame(np.zeros((new_node_A.shape[0], 4)),columns = ['Atom_label','Residue','Res_number','Note'])
        df_left['Atom_label'] = Metal_file['Atom_label']
        df_left['Residue'] = Metal_file['Residue']
        df_left['Res_number'] = Metal_count
        df_left['Note'] = Metal_file['Note']
        df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])
        df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
        df_node = pd.concat([df_node,df],ignore_index=True, join = 'outer')
        Metal_count += 1
    for i in group_B:
        new_positions=new_node_B+i
        df_left = pd.DataFrame(np.zeros((new_node_B.shape[0], 4)),columns = ['Atom_label','Residue','Res_number','Note'])
        df_left['Atom_label'] = Metal_file['Atom_label']
        df_left['Residue'] = Residue_name
        df_left['Res_number'] = Metal_count
        df_left['Note'] = Metal_file['Note']
        df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])
        df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
        df_node = pd.concat([df_node,df],ignore_index=True, join = 'outer')
        Metal_count += 1
    
    return df_node

def get_box_dimension(file):
    x1,x2,y1,y2,z1,z2 = get_box(file)
    dx,dy,dz = abs(x1-x2), abs(y1-y2), abs(z1-z2)
    dimension = [str(round(dx,5)), str(round(dy,5)),str(round(dz,5))]
    s = ' '.join(dimension)
    return(s)

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
            value_x = float(values[4])/10 #x      
            value_y = float(values[5])/10 #y
            value_z = float(values[6])/10 #z
            #value11 = values[6] #note
            # Format the values using the specified format string
            formatted_line = "%5d%-5s%5s%5d%8.4f%8.4f%8.4f" % (
                        value_resnumber, value_resname, value_label, value_atom_number, value_x, value_y, value_z) 
            newgro.append(formatted_line+'\n')        

        tail = get_box_dimension(df_all.iloc[4:6])+'\n'
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
            value_x = float(values[4]) #x      
            value_y = float(values[5]) #y
            value_z = float(values[6]) #z
            #value11 = values[6] #note
            # Format the values using the specified format string
            formatted_line = "%-5s%8.3f%8.3f%8.3f" % (
                        value_label, value_x, value_y, value_z
            )        
            newxyz.append(formatted_line+'\n')        

        fp.writelines(newxyz)


textbook_Metal_file = readpdb('Al_Al.pdb')

pAl1,pAl2,pAl3 = (textbook_Metal_file.loc[15, ['x','y','z']].to_numpy(),
                  textbook_Metal_file.loc[8, ['x','y','z']].to_numpy(),
                  textbook_Metal_file.loc[29, ['x','y','z']].to_numpy())    

pAl1_1,pAl1_2,pAl1_3 =(textbook_Metal_file.loc[16, ['x','y','z']].to_numpy(),
                          textbook_Metal_file.loc[17, ['x','y','z']].to_numpy(),
                          textbook_Metal_file.loc[18, ['x','y','z']].to_numpy()) 
   
solution_1_2,arr_1_2 = find_solution(pAl1,pAl2,pAl1_1,pAl1_2,pAl1_3)
solution_1_3,arr_1_3 = find_solution(pAl1,pAl3,pAl1_1,pAl1_2,pAl1_3)
#print(solution_1_2,solution_1_3)
#print(np.dot(solution_1_2,arr_1_2),np.dot(solution_1_3,arr_1_3))

Metal_file=readpdb('test.pdb')
axis1 = np.array([1,0,0])
axis2 =  get_axis2(solution_1_2,arr_1_2,solution_1_3,arr_1_3)

axis3 = np.cross(axis1,axis2)


point_Al = Metal_file.loc[0, ['x','y','z']].to_numpy()
p1,p2,p3 = (Metal_file.loc[1, ['x','y','z']].to_numpy()- point_Al,
                                    Metal_file.loc[2, ['x','y','z']].to_numpy()- point_Al,
                                    Metal_file.loc[3, ['x','y','z']].to_numpy()- point_Al )     
p1,p2,p3=normalize_vector(p1),normalize_vector(p2),normalize_vector(p3)                         

arr = np.vstack((p1,p2,p3))
V1,V2 = np.dot(solution_1_2,arr),np.dot(solution_1_3,arr)
V1,V2 = normalize_vector(V1),normalize_vector(V2)

Al_node = Metal_file.loc[:,['x','y','z']].to_numpy() - point_Al  #MOVE center (Al this case) to (0,0,0)
q1 = calculate_q_rotation_with_vectors(V1,axis1) #TODO:
q_V2 = quaternion.from_vector_part(V2)
new_q_V2 = q1*q_V2
new_V2 = quaternion.as_vector_part(new_q_V2)
angle = calculate_angle_rad(axis1,new_V2,axis2)
#print(new_V2,axis2,angle)
#q2 = calculate_q_rotation_with_axis_degree(axis1,angle)
q2 = calculate_q_rotation_with_vectors(new_V2,axis2)
#q3 = quaternion.from_float_array([0,0,0,-1])
#dy dz rotate pi
q3 = calculate_q_rotation_with_axis_degree(axis2,np.pi)*calculate_q_rotation_with_axis_degree(axis3,np.pi)
q_A = q2*q1
q_B = q3*q2*q1

new_node_A = get_rotate_array(Al_node,q_A)
new_node_B = get_rotate_array(Al_node,q_B)

x_num,y_num,z_num,dx_value,dy_value,dz_value = 100,5,5,3.31,16.8,15.94
group_A,group_B = group_points_AB(x_num,y_num,z_num,dx_value,dy_value,dz_value)


linker_cut_count,Residue_name = 1,'Al'
df_node = calculate_node(Metal_file,linker_cut_count,Residue_name,group_A,group_B,new_node_A,new_node_B)
df_node.to_csv('node.txt',header=None,sep='\t',index=False)

outgro(df_node,'node',0)
outxyz('node',0)
endTime = datetime.datetime.now()
print('\n'+"Time cost (s):   "+str(endTime-startTime))