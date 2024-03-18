import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def groupB_one_step_to_groupA(first_A,group_B,dx,points):
    '''same as groupA_one_step_groupB'''
    point_dx = group_B+dx
    group_A = first_A
    for i in range(group_B.shape[0]):
        point = np.reshape(point_dx[i], (1, 3))
        group_A = np.concatenate((group_A,point),axis = 0) if np.all(points == point, axis=1).any() else group_A
        group_A = np.unique(group_A,axis=0)
    return group_A

x_num,y_num,z_num,dx_value,dy_value,dz_value = 3,4,2,1,1,0.1

points=points_generator(x_num,y_num,z_num,dx_value,dy_value,dz_value )
O = np.array([[0,0,0]])
dx = dx_value*np.array([[1,0,0]]) #dx_value works as a scalar
dy = dy_value*np.array([[0,1,0]])
dz = dz_value*np.array([[0,0,1]])

first_A = O
first_B = O+dx
group_A = first_A
group_B = first_B

for i in range(points.shape[0]):
   if group_A.shape[0]+group_B.shape[0] == points.shape[0]:
      break
   else:
      group_B_dx = groupA_one_step_to_groupB(first_B,group_A,dx,points)
      group_B_dy = groupA_one_step_to_groupB(first_B,group_A,dy,points)
      group_B_dz = groupA_one_step_to_groupB(first_B,group_A,dz,points)
      group_B_dxyz = [group_B,group_B_dx,group_B_dy,group_B_dz]
      group_B = np.concatenate([arr for arr in group_B_dxyz if arr is not None],axis = 0) 
      group_B = np.unique(group_B,axis=0)
      group_A_dx = groupB_one_step_to_groupA(first_A,group_B,dx,points)
      group_A_dy = groupB_one_step_to_groupA(first_A,group_B,dy,points)
      group_A_dz = groupB_one_step_to_groupA(first_A,group_B,dz,points)
      group_A_dxyz = [group_A,group_A_dx,group_A_dy,group_A_dz]
      group_A = np.concatenate([arr for arr in group_A_dxyz if arr is not None],axis = 0)
      group_A = np.unique(group_A,axis=0)
   
print(group_A.shape,group_B.shape)

# Creating figure
fig = plt.figure()
ax = plt.axes(projection ="3d")
# Creating plot
ax.scatter3D(group_A[:,0],group_A[:,1],group_A[:,2] ,color = "green")
ax.scatter3D(group_B[:,0],group_B[:,1],group_B[:,2] ,color = "blue")
plt.title("simple 3D scatter plot")
 
# show plot
plt.show()