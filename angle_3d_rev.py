import numpy as np

def angle_val(angle):
    a=np.array(angle)
    v1=np.subtract(a[1],a[0])
    v2=np.subtract(a[2],a[1])
    return np.acos(np.matmul(v1,v2.T)/(np.matmul(v1,v1.T)*np.matmul(v2,v2.T))**0.5)

def angle_derivative_vector(angle,vector):
    a=np.array(angle)
    v=np.array(vector)

    derivative_len=0.00001

    v_all=np.concatenate(v)
    v_len=np.linalg.norm(v_all)
    v_der=np.multiply(derivative_len/v_len,v)
    return (angle_val(np.subtract(a,v_der)) - angle_val(np.add(a,v_der)))/(2*derivative_len)

def find_vectors(angle):
    a=np.array(angle)

    v_ac=a[2]-a[0]
    v_ab=a[1]-a[0]
    v_n=np.cross(v_ac,v_ab)

    a_sqr=np.matmul(a[0],a[0].T)
    b_sqr=np.matmul(a[1],a[1].T)
    c_sqr=np.matmul(a[2],a[2].T)

    print(v_n)
 

    o=np.linalg.solve(np.array([v_ac,v_ab,v_n]),np.array([(c_sqr-a_sqr)*0.5, (b_sqr-a_sqr)*0.5, np.matmul(a[0],v_n.T)]))
    
    v_ob=np.subtract(a[1],o)
    v_a_n=np.cross(v_ab,v_n)
    v_c_n=np.cross(np.subtract(a[1],a[2]),v_n)

    v0=[0,0,0]

    deriv_0=angle_derivative_vector(a,[v_a_n,v0,v0])
    deriv_1=angle_derivative_vector(a,[v0, v_ob,v0])
    deriv_2=angle_derivative_vector(a,[v0,v0,v_c_n])
    v_deriv_0=np.multiply(deriv_0/np.linalg.norm(v_a_n),v_a_n)
    v_deriv_1=np.multiply(deriv_1/np.linalg.norm(v_a_n),v_ob)
    v_deriv_2=np.multiply(deriv_2/np.linalg.norm(v_a_n),v_c_n)

    return v_deriv_0,v_deriv_1,v_deriv_2

if __name__ == "__main__":
    a=[2.154,0,0]
    b=[0,0,0]
    c=[0,0,2.154]
    print(np.array(find_vectors([a,b,c])))