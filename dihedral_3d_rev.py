import numpy as np

def dihedral_val(dihedral):
    d=np.array(dihedral)
    u1=np.subtract(d[1],d[0])
    u2=np.subtract(d[2],d[1])
    u3=np.subtract(d[3],d[2])
    u1pu2=np.cross(u1,u2)
    u2pu3=np.cross(u2,u3)
    arg1=np.matmul(u2,np.cross(u1pu2,u2pu3).T)
    arg2=np.linalg.norm(u2)*np.matmul(u1pu2,u2pu3.T)
    return np.arctan2(arg1, arg2)

def dihedral_derivative_vector(dihedral,vector):
    d=np.array(dihedral)
    v=np.array(vector)

    derivative_len=0.00001

    v_all=np.concatenate(v)
    v_len=np.linalg.norm(v_all)
    v_der=np.multiply(derivative_len/v_len,v)
    #добавить что-то про одолистность вблизи -pi и pi 
    return (dihedral_val(np.subtract(d,v_der)) - dihedral_val(np.add(d,v_der)))/(2*derivative_len)

def find_vectors(dihedral):
    d=np.array(dihedral)

    v_ab=np.subtract(d[1],d[0])
    v_cd=np.subtract(d[3],d[2])
    v_n=np.subtract(d[2],d[1])#v_bc

    v0=[0,0,0]

    v_a_n=np.cross(v_ab,v_n)
    v_d_n=np.cross(v_cd,v_n)
    
    deriv_0=dihedral_derivative_vector(d,[v_a_n,v0,v0,v0])
    deriv_3=dihedral_derivative_vector(d,[v0,v0,v0,v_d_n])
    v_deriv_0=np.multiply(deriv_0/np.linalg.norm(v_a_n), v_a_n)
    v_deriv_3=np.multiply(deriv_3/np.linalg.norm(v_d_n), v_d_n)

    vx,vy,vz=[1,0,0],[0,1,0],[0,0,1]

    v_deriv_1=np.array([dihedral_derivative_vector(d,[v0,vx,v0,v0]), dihedral_derivative_vector(d,[v0,vy,v0,v0]), dihedral_derivative_vector(d,[v0,vz,v0,v0])])
    v_deriv_2=np.array([dihedral_derivative_vector(d,[v0,v0,vx,v0]), dihedral_derivative_vector(d,[v0,v0,vy,v0]), dihedral_derivative_vector(d,[v0,v0,vz,v0])])

    return v_deriv_0,v_deriv_1,v_deriv_2,v_deriv_3

if __name__=="__main__":
    a=[1,0,0]
    b=[0,0,0]
    c=[0,1,0]
    d=[0,1,1]
    print(np.array(find_vectors([a,b,c,d])))