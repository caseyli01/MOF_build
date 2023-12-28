import numpy as np
import pandas as pd
import re

def get_box(file):
    x, y, z = [], [], []
    x, y, z = file['x'],file['y'],file['z']
    x1,y1,z1 = min(x), min(y), min(z)
    x2,y2,z2 = max(x), max(y), max(z)
    return(x1,x2,y1,y2,z1,z2)

def if_points_out_of_box(df,frame):
    x1,x2,y1,y2,z1,z2 = get_box(frame)[0],get_box(frame)[1],get_box(frame)[2],get_box(frame)[3],get_box(frame)[4],get_box(frame)[5]
    results = []
    for i in range(df.shape[0]):
        A = df.iloc[i,:]
        result_X = df.iloc[i,0] < x1 or df.iloc[i,0] > x2
        result_Y = df.iloc[i,1] < y1 or df.iloc[i,1] > y2
        result_Z = df.iloc[i,2] < z1 or df.iloc[i,2] > z2
        result = result_X or result_Y or result_Z
        results.append(result)
    return any(results)


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
    #newpath = os.path.abspath ( '')+'/'+str(outputfile)+'/'    # input file
    #os.makedirs(newpath,exist_ok=True) 
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
    

    data = pd.read_csv(outputfile+'.txt',delim_whitespace=True,names=['Atom_label','Residue','Res_number','x','y','z','Note'])
    return data


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

def are_points_collinear_3d(pointA, pointB, pointC):
    a1,a2,a3 = pointA
    b1,b2,b3 = pointB
    c1,c2,c3 = pointC
    x1, y1, z1 = (b1-a1),(b2-a2),(b3-a3)
    x2, y2, z2 = (c1-a1),(c2-a2),(c3-a3)
    #print("linear")
    return x2*y1== y2*x1 and y2*z1==y1*z2

def calculate_normal_vector(pointA, pointB, pointC):
    vectorAB = pointB - pointA
    #vectorAB=vectorAB.round(3)
    vectorAC = pointC - pointA
    #vectorAC=vectorAC.round(3)
    normal_vector = np.cross(vectorAB, vectorAC)
    unit_normal_vector = normal_vector/np.linalg.norm(normal_vector)
    return unit_normal_vector

def calculate_rotation_angle(point_A, point_B,point_C):
    vector_AB = point_B - point_A
    vector_AC = point_C - point_A
    # Define vectors A and B
    dot_AB_AC= np.dot(vector_AB, vector_AC)
    # Calculate the magnitudes of A and B
    magnitude_AB = np.linalg.norm(vector_AB)
    magnitude_AC = np.linalg.norm(vector_AC)
    
    cos_theta = dot_AB_AC/(magnitude_AB*magnitude_AC)
    theta_rad = np.arccos(cos_theta)
    #theta_degrees = np.degrees(theta)
    return theta_rad

def calculate_rotation_matrix(theta,rotation_axis):# Define rotation angle and axis
   
    # Calculate the skew-symmetric cross-product matrix K
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    
    # Calculate the rotation matrix using Rodrigues' formula
    R = np.identity(3) + np.sin(theta) * K + (1 - np.cos(theta)) * np.dot(K, K)
    R=np.round(R,2)
    return R


def find_points_in_radius(points, center, radius):
    points_identity = points.iloc[:, [0,1,2,5]]
    points_positions = np.array(points.loc[:, ['x','y','z']])
    center = np.array(center)
    distances = np.linalg.norm(points_positions - center, axis=1)
    indices = np.where(distances <= radius)[0]
    #result = pd.concat([points_identity.iloc[indices,0], points[indices]], axis=1)
    return points[indices]

def find_points_in_cutoff(points, center, side_atom, radius):
    points_positions = np.array(points.loc[:, ['x','y','z']])
    center = np.array(center)
    side_atom =np.array(side_atom)
    distances = np.linalg.norm(points_positions - center, axis=1)
    indices = np.where(distances <= radius)[0]
    out_points_center = points.iloc[indices]
    #out_points_center.loc[:, ['x','y','z']] = 2*center-out_points_center.loc[:, ['x','y','z']] #over A to outchain
    out_points = out_points_center
    out_points.iloc[:,1] = 'CUT'
    return out_points


def find_farthest_point_in_cutoff(points, center):
    points_positions = np.array(points.loc[:, ['x','y','z']])
    center = np.array(center)
    distances = np.linalg.norm(points_positions - center, axis=1)
    farthest_atom_index = np.argmax(distances)
    farthest_atom = points.iloc[farthest_atom_index]
    #print(farthest_atom)
    return farthest_atom


def calculate_user_deifined_termination(He,Cut_in_linker,input,index,d):
    He_points = np.array(He.loc[:, ['x','y','z']])
    Y = He_points[0]
    farthest_atom = find_farthest_point_in_cutoff(Cut_in_linker,Y)
    FA = np.array(farthest_atom.loc[['x','y','z']])
    Y_FA= FA-Y
    oldC = np.array(input.loc[index-1, ['x','y','z']])
    newC= FA+d*normalize_vector(Y_FA)
    translation = newC-oldC
    extra = input.loc[:, ['x','y','z']]+translation
    input.loc[:,['x','y','z']] = extra
    input.loc[:,'Residue'] = 'CUT'
    #print(input)
    return input



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
                #vAB =  [float(str(val)[:str(val).index('.') + 3]) for val in v]
                #set A to (0,0,0)
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
            #vAB = np.trunc(vAB* 100) / 100
            #set A to (0,0,0)
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
        #print(filtered_vector,filtered_vector.index)
        term_vector.append(filtered_vector.to_numpy())
        true_indices.append(filtered_vector.index)
        #print([filtered_vector,true_indices])
    
    return [term_vector,true_indices]


def unique_rotation_matrix(unique_vector,point_A,point_B):
    #set A to (0,0,0)
    rotation_matrixes = []
    for i in range(len(unique_vector)):
        point_A_0 = point_A -point_A
        AB_S1 = point_B -point_A #in linker file
        AB_S2 = unique_vector.iloc[i,:] # in frame file
        
        if are_points_collinear_3d(point_A_0,AB_S1,AB_S2)!=True:
            rotation_axis = calculate_normal_vector(point_A_0,AB_S1,AB_S2) 
            rotation_angle = -1*calculate_rotation_angle(point_A_0,AB_S1,AB_S2)
            rotation_matrix = calculate_rotation_matrix(rotation_angle,rotation_axis)
            rotation_matrixes.append(rotation_matrix)
        else:
            rotation_matrix = np.array([[1,0,0],[0,1,0],[0,0,1]])
        
    return rotation_matrixes

def get_zoomedframe(points_n,frame,zoom):
    df1 = pd.DataFrame(columns = ['Atom_label','Residue','Res_number','x','y','z','Note'])
    for i in range(points_n):
        df1.loc[i,'Atom_label'] = frame.loc[i,'Atom_label']
        df1.loc[i,'Residue'] = frame.loc[i,'Residue']
        df1.loc[i,'Res_number'] = str(i+1)
        df1.loc[i,'Note'] = frame.loc[i,'Note']
        df1.loc[i,'x'] = frame.loc[i,'x'] *zoom
        df1.loc[i,'y'] = frame.loc[i,'y'] *zoom
        df1.loc[i,'z'] = frame.loc[i,'z'] *zoom
    return df1


def calculate_Zr_octahedral(point_A,points_n,MM_l,df1,Metal_file,linker_cut_count,Residue_name):
#octahedral
    octa_vertices = np.array([
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.]
    ])
    Metal_count = linker_cut_count
    df_octa = pd.DataFrame()
    point_Zr = Metal_file.loc[0, ['x','y','z']]
    d_Zr_center = np.linalg.norm(point_Zr-point_A)
    for i in range(points_n-1):
        for j in range(i,points_n):
            x = round(length_square(df1.loc[i],df1.loc[j]))
            if ( x == round(MM_l)):
                point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
                point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
                #set A to (0,0,0)
                point_A_0 = np.array([0.,0.,0.])
                AB_S1 = 1/2*(octa_vertices[0]+octa_vertices[2])
                #AB_S1 = octa_vertices[0]
                AB_S2 = point_B_frame -point_A_frame
                if are_points_collinear_3d(point_A_0,AB_S1,AB_S2)!=True:
                    rotation_axis = calculate_normal_vector(point_A_0,AB_S1,AB_S2)
                    rotation_angle = -1*calculate_rotation_angle(point_A_0,AB_S1,AB_S2)
                    rotation_matrix = calculate_rotation_matrix(rotation_angle,rotation_axis)
                    #print(i,j,x,AB_S1,AB_S2)
                    #print("matrix")
                    #print(rotation_matrix)
                    new_positions=(d_Zr_center*np.dot(octa_vertices, rotation_matrix))+point_A_frame
                    #print("new position")
                    #print(i,j)
                    break
                else: 
                    new_positions=(d_Zr_center*octa_vertices)+point_A_frame
                    #print("no--------matrix")
                    #print(i,j,x,AB_S1,AB_S2)
                    break
                   
        df_left = pd.DataFrame(np.zeros((6, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
        
        df_left['Atom_label'] = Metal_file.loc[0,'Atom_label']
        df_left['Residue'] = Residue_name
        df_left['Res_number'] = Metal_count
        df_left['Note'] = Metal_file.loc[0,'Note']
        #print(i,j,Metal_count)

        df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])
        df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
        df_octa = pd.concat([df_octa,df],ignore_index=True, join = 'outer')
        #print(df_left,df_right,df_mof)
        Metal_count += 1
                #print(i,j)
    #print(df_left,df_right,new_positions)
    for j in range(points_n):
            i = points_n-1
            x = round(length_square(df1.loc[i],df1.loc[j]))
            if ( x == round(MM_l)):
                point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
                point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
                #set A to (0,0,0)
                point_A_0 = np.array([0.,0.,0.])
                #AB_S1 = octa_vertices[0]
                AB_S1 = 1/2*(octa_vertices[0]+octa_vertices[2])
                AB_S2 = point_B_frame -point_A_frame
                if are_points_collinear_3d(point_A_0,AB_S1,AB_S2)!=True:
                    rotation_axis = calculate_normal_vector(point_A_0,AB_S1,AB_S2)
                    rotation_angle = -1*calculate_rotation_angle(point_A_0,AB_S1,AB_S2)
                    rotation_matrix = calculate_rotation_matrix(rotation_angle,rotation_axis)
                    #print(i,j,x,AB_S1,AB_S2)
                    #print("matrix")
                    #print(rotation_matrix)
                    new_positions=(d_Zr_center*np.dot(octa_vertices, rotation_matrix))+point_A_frame
                    #print("new position")
                    #print(new_positions)
                    break
                else: 
                    new_positions=(d_Zr_center*octa_vertices)+point_A_frame
                    #print(i,j,x,AB_S1,AB_S2)
                    #print("no--------matrix")
                    break

    df_left = pd.DataFrame(np.zeros((6, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
        
    df_left['Atom_label'] = Metal_file.loc[0,'Atom_label']
    df_left['Residue'] = Residue_name
    df_left['Res_number'] = Metal_count
    df_left['Note'] = Metal_file.loc[0,'Note']
    #print(new_positions)

    df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])
    df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
    df_octa = pd.concat([df_octa,df],ignore_index=True, join = 'outer')
    #print(df_left,df_right,df_mof)
    Metal_count += 1


    return df_octa


def calculate_angle(vector_AB, vector_AC):
    # Define vectors A and B
    dot_AB_AC= np.dot(vector_AB, vector_AC)
    # Calculate the magnitudes of A and B
    magnitude_AB = np.linalg.norm(vector_AB)
    magnitude_AC = np.linalg.norm(vector_AC)
    cos_theta = dot_AB_AC/(magnitude_AB*magnitude_AC)
    theta_rad = np.arccos(cos_theta)
    theta_degrees = np.degrees(theta_rad)

    return theta_degrees


def find_correspond_octa_tops(HE_HE,octa):
    Zr_octa=[]
    for j in octa:
        angle = calculate_angle(HE_HE,j)
        if 44<angle<46:
            Zr_octa.append(j)
    return Zr_octa


def get_cubeface(array1):
    length = []
    for i in range(1):
        for j in  range(len(array1)):
            v= array1[i]-array1[j]
            length.append(round(np.linalg.norm(v),1))

    length = list(set(length))
    length.sort()

    b = []
    c = []
    for i in range(1):
        for j in  range(len(array1)):
            v= array1[i]-array1[j]
            if 1.1*length[1]>np.linalg.norm(v) > 0.99*length[1]:
                #b.append(i)
                b.append(j)
            if 1.1*length[3]>np.linalg.norm(v) > 0.99*length[3]:
                #b.append(i)
                c.append(j)
    topface=array1[b]
    bottomface=array1[c]
    return(topface,bottomface)

def find_cube(pointA,bottomface):
    diagonal_length = []
    bottomindex = []
    for i in range(len(bottomface)):
        v= bottomface[i] - pointA
        l=np.linalg.norm(v)
        diagonal_length.append(l)
        bottomindex.append(i)
    Gindex = diagonal_length.index(max(diagonal_length))  #body diagonal
    Eindex = diagonal_length.index(min(diagonal_length))   #vertical
    pointG=bottomface[Gindex]
    pointE=bottomface[Eindex]
    HFindex = [x for x in bottomindex if x not in [Gindex,Eindex]] # face diagonal
    pointH = bottomface[HFindex[0]]
    pointF = bottomface[HFindex[1]]
    return(pointG,pointE,pointH,pointF)

def normalize_vector(v):
    norm_v=v/np.linalg.norm(v)
    return norm_v

def get_metal_octahedral(array1):

    topface,bottomface = get_cubeface(array1)
    pointA=topface[0]
    pointG,pointE,pointH,pointF = find_cube(pointA,bottomface)
    pointA,pointC,pointB,pointD = find_cube(pointG,topface)

    #AC AH AF GD GB GE
    AC = normalize_vector(pointA+pointC)
    AH = normalize_vector(pointA+pointH)
    AF = normalize_vector(pointA+pointF)
    GD = normalize_vector(pointG+pointD)
    GB = normalize_vector(pointG+pointB)
    GE = normalize_vector(pointG+pointE)

    octahedral = np.array([AC, AH, AF, GD, GB,GE])
    #print(octahedral.shape)
    return(octahedral)


def find_an_diagnol_inframe(df1,MM_l):
    i = 1
    j = 2
    edge = []
    while (j < df1.shape[0] and len(edge) <3):
        x = round(length_square(df1.loc[i],df1.loc[j]))
        if ( x == round(MM_l)):
            point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
            point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
            AB_S2 = point_B_frame - point_A_frame
            edge.append(AB_S2)
            j +=1
   
    return edge[0]+edge[1]

def find_points_in_octa(octa):
    lens = []
    index = []
    for i in range(len(octa)):
        for j in range(len(octa)):
            v= octa[j]-octa[i]
            l=np.linalg.norm(v)
            l=round(l,2)
            lens.append(l)
            index.append([i,j])
    ls=list(set(lens))
    ls.sort()
    EFindex = lens.index(ls[-1])  #body diagonal
    ACindex,BDindex = lens.index(ls[-2])   #vertical
    E=index[EFindex][0]
    F=index[EFindex][1]
    A=index[ACindex][0]
    B=index[BDindex][0]
    #rotation_axis=EF, OA-->OB

def calculate_cos(vector1,vector2):
        dot_product = np.dot(vector1, vector2)
        magnitude_vector1 = np.linalg.norm(vector1)
        magnitude_vector2 = np.linalg.norm(vector2)
        cosine_angle = dot_product / (magnitude_vector1 * magnitude_vector2)
        angle_rad = np.arccos(np.clip(cosine_angle, -1.0, 1.0))    
        angle_deg = np.degrees(angle_rad)
        return angle_deg

def rotated_node(octa,array1):
# find 2 vector who doesn't parallel to any vector in unique vectors
    angle = []
    EF = []
    EFindex = []
    for i in range(len(octa)):
        angle = []
        for j in range(len(array1)):
            cos = calculate_cos(octa[i],array1[j])
            angle.append(cos)
        if 0 not in angle:
            EF.append(octa[i])
            EFindex.append(i)
    E=EF[0]
    F=EF[1]

    ABCDface=np.delete(octa, EFindex, axis=0)

    #find ABCD
    A=ABCDface[0]
    BD=[]
    for j in range(1,len(ABCDface)):
        cos = calculate_cos(A,ABCDface[j])
        if 179<cos<181:
            C=ABCDface[j]
        else:
            BD.append(ABCDface[j])
    B=BD[0]
    D=BD[1]
    rotation_axis=E-F
    A1B=0.707*(A+B)
    A1D=0.707*(A+D)
    C1B=0.707*(C+B)
    C1D=0.707*(C+D)

    newocta=np.array([E,F,A1B,A1D,C1B,C1D])
    
    return (newocta)



def calculate_MOF_linker(point_A,point_B,points_n,MM_l,df1,L_filename,Zr_linker,count):
    df_mof = df1 
    #TODO:df_mof = pd.DataFrame() but need to be sure the columns positions
    residue_count = count
    unique_vector = np.array(search_unique_vector(df1,MM_l))
    octa_vertices = get_metal_octahedral(unique_vector)
    octa=rotated_node(octa_vertices,unique_vector)
    for i in range(points_n):
#for i in range(6,7):
        for j in range(i,points_n):
            x = round(length_square(df1.loc[i],df1.loc[j]))
            if ( x == round(MM_l)):
                point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
                point_B_frame =  np.asarray(df1.loc[j,['x','y','z']],dtype = float)
                #set A to (0,0,0)
                point_A_0 = point_A -point_A
                AB_S1 = point_B -point_A
                AB_S2 = point_B_frame -point_A_frame

                rotation_axis = calculate_normal_vector(point_A_0,AB_S1,AB_S2)
                rotation_angle = -1*calculate_rotation_angle(point_A_0,AB_S1,AB_S2)
                rotation_matrix = calculate_rotation_matrix(rotation_angle,rotation_axis)
                linker_positions=np.dot(L_filename.loc[:,['x','y','z']].values-point_A, rotation_matrix)
                # use Zr_linker to locate self rotation matrix 
                
                newZr_inlinker = np.dot(Zr_linker.loc[:,['x','y','z']].values-point_A, rotation_matrix)               
                vZr_inlinker = newZr_inlinker[0]-newZr_inlinker[1]
                topZr_inocta = find_correspond_octa_tops(AB_S2,octa)
                #print(topZr_inocta)
                vZr_incluster = topZr_inocta[0]-topZr_inocta[1]
                self_rotation_axis = calculate_normal_vector(point_A_0,vZr_incluster,vZr_inlinker)
                self_rotation_angle =  calculate_rotation_angle(point_A_0,vZr_incluster,vZr_inlinker)
                #print(self_rotation_angle,self_rotation_axis)        
                self_rotation_matrix = calculate_rotation_matrix(self_rotation_angle,self_rotation_axis) 
                #print(self_rotation_matrix)
                #new_positions = np.dot(linker_positions-newZr_inlinker[1],self_rotation_matrix)+newZr_inlinker[1]+point_A_frame                
                midZr_octa= 0.5*(topZr_inocta[0]+newZr_inlinker[1])
                new_positions = np.dot(linker_positions-midZr_octa,self_rotation_matrix)+midZr_octa+point_A_frame


                df_left = pd.DataFrame(columns = ['Atom_label','Residue','Res_number','Note'])
                df_left.loc[:,'Atom_label'] = L_filename.loc[:,'Atom_label']
                df_left.loc[:,'Residue'] = L_filename.loc[:,'Residue']
                df_left.loc[:,'Res_number'] = residue_count
                df_left.loc[:,'Note'] = L_filename.loc[:,'Note']
                
                df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])
                df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
                df_mof = pd.concat([df_mof,df], ignore_index=True,keys=['df_mof', 'df'], join = 'outer')
                #print(df_left,df_right,df_mof)
                residue_count += 1
    return df_mof


def calculate_MOF_CUT_TERM(point_A,point_B,points_n,MM_l,df1,PMMP_file,Metal_file,Zr_linker,count,Residue_name):
    unique_vector = np.array(search_unique_vector(df1,MM_l))
    octa_vertices = get_metal_octahedral(unique_vector)

    point_Zr = Metal_file.loc[0, ['x','y','z']]
   
    d_Zr_center = np.linalg.norm(point_Zr-point_A) 
    #for i in range(6,7):
    octa=rotated_node(octa_vertices,unique_vector)
    zoomed_octa=d_Zr_center*octa

    unique_vector = search_unique_vector(df1,MM_l)
    term_vector = filtered_term_vector(df1,MM_l)[0]
    true_indices = filtered_term_vector(df1,MM_l)[1]
    rotation_matrixes = unique_rotation_matrix(unique_vector,point_A,point_B)
    residue_Term_count = count
    df_outpoints = pd.DataFrame()
    for i in range(points_n):
        point_term_vector = term_vector[i]

        for j in range(len(true_indices[i])):
                rotation_matrix_index = true_indices[i]
                point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
                #set A to (0,0,0)
                point_A_0 = point_A -point_A
                AB_S1 = point_B -point_A #in linker file
                AB_S2 = unique_vector.iloc[j,:] # in frame file
                if are_points_collinear_3d(point_A_0,AB_S1,AB_S2)!=True:
                    rotation_matrix = rotation_matrixes[rotation_matrix_index[j]]
                    linker_positions=np.dot(PMMP_file.loc[:, ['x','y','z']].values-point_A, rotation_matrix)
                    newZr_inlinker = np.dot(Zr_linker.loc[:,['x','y','z']].values-point_A, rotation_matrix)
                else: 
                    #linker can tranlate to edge without rotation
                    linker_positions=PMMP_file.loc[:, ['x','y','z']].values-point_A
                    newZr_inlinker = Zr_linker.loc[:,['x','y','z']].values-point_A
              
                           
                vZr_inlinker = newZr_inlinker[0]-newZr_inlinker[1]
                topZr_inocta = find_correspond_octa_tops(point_term_vector[j],zoomed_octa)
                #print(topZr_inocta)
                vZr_incluster = topZr_inocta[0]-topZr_inocta[1]
                self_rotation_axis = calculate_normal_vector(point_A_0,vZr_incluster,vZr_inlinker)
                self_rotation_angle =  calculate_rotation_angle(point_A_0,vZr_incluster,vZr_inlinker)    
                self_rotation_matrix = calculate_rotation_matrix(self_rotation_angle,self_rotation_axis) 
                #new_positions = np.dot(linker_positions-newZr_inlinker[1],self_rotation_matrix)+newZr_inlinker[1]+point_A_frame                
                new_positions = np.dot(linker_positions-newZr_inlinker[1],self_rotation_matrix)+topZr_inocta[1]+point_A_frame
               

                df_right = pd.DataFrame(new_positions,columns = ['x','y','z'])       
                df_left = pd.DataFrame(columns = ['Atom_label','Residue','Res_number','Note'])
                df_left.loc[:,'Atom_label'] = PMMP_file.loc[:,'Atom_label'] 
                df_left.loc[:,'Residue'] = Residue_name
                #df_left.loc[:,'Residue'] = PMMP_file.loc[:,'Residue']
                df_left.loc[:,'Res_number'] = residue_Term_count
                df_left.loc[:,'Note'] = PMMP_file.loc[:,'Note']
                #df_left = df_left.reset_index(drop = True)
                df = pd.concat([df_left,df_right],axis = 1) 
                residue_Term_count += 1
                df_outpoints = pd.concat([df_outpoints,df],ignore_index=True, join = 'outer')
    return df_outpoints



def calculate_dummy_D1(He,Zr1,Zr2):
    Y = 0.5*(Zr1+Zr2)
    HeZr1 = Zr1-He
    HeY = Y-He
    YZr1=Zr1-Y
    x= normalize_vector(HeY)
    y= normalize_vector(YZr1)
    
    O1 = 4.0*x+1.1*y+He ##calculated from hMOF-1003162 structure

           
    #Zr1=(A1,B1); O(A2,B2); dummyD1(A3,B3),opposite dummyD2(A4,B4) related to He as HeZr1,HeO1,HeD1
    #A1 = np.dot(HeZr1,x.reshape(3,1))[0]
    #B1 = np.dot(HeZr1,y.reshape(3,1))[0]

    #HeO1 = O1-He
    #A2 = np.dot(HeO1,x.reshape(3,1))[0]
    #B2 = np.dot(HeO1,y.reshape(3,1))[0]

    Zr1O1=O1-Zr1
    ###this is angstrom unit so it is 1A
    #A3 = A1+normalize_vector(Zr1O1)
    #B3 = B1+normalize_vector(Zr1O1)

    #A4 = A1-normalize_vector(Zr1O1)
    #B4 = B1-normalize_vector(Zr1O1)
    

    HeD1 = HeZr1+normalize_vector(Zr1O1)
    #HeD2 = HeZr1+normalize_vector(Zr1O2)

    return(HeD1)

    

def get_dummy_atoms(He,Zr_nodes,Zr1):
    dummy = []
    for j in range(len(Zr_nodes)):
        if 5<calculate_angle(Zr1,Zr_nodes[j]) <175:
            D1 = calculate_dummy_D1(He,Zr1,Zr_nodes[j])
            dummy.append(D1)
    #calculate D2
    Zr2_pairs = []
    for j in range(len(Zr_nodes)):
        if 5<calculate_angle(Zr1,Zr_nodes[j]) <175:
            Zr2_pairs.append(Zr_nodes[j])
    vZr1_Zr2 = Zr2_pairs -Zr1
    i_index, j_index = [], []
    for i in range(len(vZr1_Zr2)):
        for j in range(i,len(vZr1_Zr2)):
            if 58<calculate_angle(vZr1_Zr2[i],vZr1_Zr2[j]) <62:
                i_index.append(i)
                j_index.append(j)
    #print(i_index,j_index)
    O2 = []
    for k in range(len(i_index)):
        O2.append(0.333*(Zr1+Zr2_pairs[i_index[k]]+Zr2_pairs[j_index[k]]))

    Zr1O2=O2-Zr1
    HeZr1 = Zr1-He
    for i in range(len(i_index)):
        dummy.append(HeZr1+normalize_vector(Zr1O2[i]))
    return dummy


def get_OOH_in_node(octa):
    E,F,A,D,B,C = octa[0],octa[1],octa[2],octa[3],octa[4],octa[5]
    O_m2 = []
    O_m2.append(0.333*(F+C+D))
    O_m2.append(0.333*(F+A+B))
    O_m2.append(0.333*(E+A+D))
    O_m2.append(0.333*(E+C+B))
    O_oh = []
    O_oh.append(0.333*(F+A+D))
    O_oh.append(0.333*(F+B+C))
    O_oh.append(0.333*(E+C+D))
    O_oh.append(0.333*(E+A+B))
    H_oh= []
    H_oh.append(O_oh[0]-0.981*calculate_normal_vector(F,A,D))
    H_oh.append(O_oh[1]+0.981*calculate_normal_vector(F,B,C))
    H_oh.append(O_oh[2]-0.981*calculate_normal_vector(E,C,D))
    H_oh.append(O_oh[3]-0.981*calculate_normal_vector(E,A,B))
    return O_m2,O_oh,H_oh



def get_framenode_octahedral(point_A,points_n,MM_l,df1,Metal_file,linker_cut_count,Residue_name,residue_path):
#octahedral
    unique_vector = np.array(search_unique_vector(df1,MM_l))
    octa_vertices = get_metal_octahedral(unique_vector)
    #print(unique_vector)
    #print(octa_vertices)
    Metal_count = linker_cut_count
    df_octa = pd.DataFrame()
    point_Zr = Metal_file.loc[0, ['x','y','z']]
   
    d_Zr_center = np.linalg.norm(point_Zr-point_A) 
    #for i in range(6,7):
    octa=rotated_node(octa_vertices,unique_vector)
    zoomed_octa=d_Zr_center*octa
    #print(zoomed_octa)
    for i in range(points_n):
        point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
        for j in range(len(zoomed_octa)):
            new_positions=zoomed_octa[j]+point_A_frame
            df_left = pd.DataFrame(np.zeros((1, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
            df_left['Atom_label'] = Metal_file.loc[0,'Atom_label']
            df_left['Residue'] = Residue_name
            df_left['Res_number'] = Metal_count
            df_left['Note'] = Metal_file.loc[0,'Note']
            #print(i,j,Metal_count)
            df_right = pd.DataFrame(new_positions.reshape(1,3),columns = ['x','y','z'])
            df = pd.concat([df_left,df_right],axis = 1, join = 'outer') 
            
            He = point_A_frame-point_A_frame
            dummy_array=get_dummy_atoms(He,zoomed_octa,zoomed_octa[j])+point_A_frame

            Ddf_left = pd.DataFrame(np.zeros((8, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
            Ddf_left['Atom_label'] = ['D1','D2','D3','D4','D5','D6','D7','D8']
            Ddf_left['Residue'] = Residue_name
            Ddf_left['Res_number'] = Metal_count
            Ddf_left['Note'] = Metal_file.loc[0,'Note']
            #print(i,j,Metal_count)
            Ddf_right = pd.DataFrame(dummy_array,columns = ['x','y','z'])
            Ddf = pd.concat([Ddf_left,Ddf_right],axis = 1, join = 'outer') 
            df_octa = pd.concat([df_octa,df,Ddf],ignore_index=True, join = 'outer')
           
            #print(df_left,df_right,df_mof)
            Metal_count += 1
                #print(i,j)
            
    ZR_residue = pd.concat([df,Ddf],ignore_index=True, join = 'outer')
    ZR_residue.to_csv(residue_path+'ZR1.txt', sep='\t', header = None)
    O_m2,O_oh,H_oh = get_OOH_in_node(zoomed_octa)
    for i in range(points_n):
        point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
        for j in range(len(O_m2)):
            new_positions=O_m2[j]+point_A_frame
            O_m2_df_left = pd.DataFrame(np.zeros((1, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
            O_m2_df_left['Atom_label'] = 'O2'
            O_m2_df_left['Residue'] = 'OM2'
            O_m2_df_left['Res_number'] = Metal_count
            O_m2_df_left['Note'] = 'O2-'
            #print(i,j,Metal_count)
            O_m2_df_right = pd.DataFrame(new_positions.reshape(1,3),columns = ['x','y','z'])
            O_m2_df = pd.concat([O_m2_df_left,O_m2_df_right],axis = 1, join = 'outer') 
            df_octa = pd.concat([df_octa,O_m2_df],ignore_index=True, join = 'outer')

            Metal_count += 1
    O_m2_df.to_csv(residue_path+'O_m2.txt', sep='\t', header = None)
    for i in range(points_n):
        point_A_frame =  np.asarray(df1.loc[i,['x','y','z']],dtype = float)
        for j in range(len(O_oh)):
            new_positions=O_oh[j]+point_A_frame
            O_oh_df_left = pd.DataFrame(np.zeros((1, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
            O_oh_df_left['Atom_label'] = 'O1'
            O_oh_df_left['Residue'] = 'OH'
            O_oh_df_left['Res_number'] = Metal_count
            O_oh_df_left['Note'] = 'O_OH'
            #print(i,j,Metal_count)
            O_oh_df_right = pd.DataFrame(new_positions.reshape(1,3),columns = ['x','y','z'])
            O_oh_df = pd.concat([O_oh_df_left,O_oh_df_right],axis = 1, join = 'outer') 

            new_positions=H_oh[j]+point_A_frame
            H_oh_df_left = pd.DataFrame(np.zeros((1, 4)),columns = ['Atom_label','Residue','Res_number','Note'])
            H_oh_df_left['Atom_label'] = 'H'
            H_oh_df_left['Residue'] = 'OH'
            H_oh_df_left['Res_number'] = Metal_count
            H_oh_df_left['Note'] = 'H_OH'
            #print(i,j,Metal_count)
            H_oh_df_right = pd.DataFrame(new_positions.reshape(1,3),columns = ['x','y','z'])
            H_oh_df = pd.concat([H_oh_df_left,H_oh_df_right],axis = 1, join = 'outer') 
            df_octa = pd.concat([df_octa,O_oh_df,H_oh_df],ignore_index=True, join = 'outer')
            Metal_count += 1
    pd.concat([O_oh_df,H_oh_df],ignore_index=True, join = 'outer').to_csv(residue_path+'OH.txt', sep='\t', header = None)

    return df_octa
