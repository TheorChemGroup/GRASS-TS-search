def mirror_fn(grad,#list 3*N, gradient
              xyzs,#list 3*N, coord i.e.[[0, 1, 1],...]
              search_DoFs,#list of reaction dofs: type, athoms, value. i.e. [["b", 1, 2, -1"],...] for stratching (-1) bond (b) between 1 and 2 athoms
              be_verbose
              ):
    #print(search_DoFs)
    import numpy as np
    from angle_3d_rev import find_vectors as f_v_a
    from dihedral_3d_rev import find_vectors as f_v_d
    m_grad=np.array(grad)
    nAthoms=len(xyzs)
    
    v0=[0.,0.,0.]
    mirror_vec=[]
    for i in range(nAthoms):
        mirror_vec.append(v0)
    mirror_vec=np.array(mirror_vec)
    
    for dof in search_DoFs:
        if dof[0]=='b':
            direction=np.subtract(xyzs[dof[1]-1], xyzs[dof[2]-1])
            direction=np.multiply(dof[3]/np.linalg.norm(direction),direction)
            mirror_vec[dof[1]-1]-=direction
            mirror_vec[dof[2]-1]+=direction
        
        if dof[0]=='a':
            vectors=f_v_a([xyzs[dof[1]-1],xyzs[dof[2]-1],xyzs[dof[3]-1]])
            mirror_vec[dof[1]-1]+=dof[4]*vectors[0]
            mirror_vec[dof[2]-1]+=dof[4]*vectors[1]
            mirror_vec[dof[3]-1]+=dof[4]*vectors[2]

        if dof[0]=='d':
            vectors=f_v_d([xyzs[dof[1]-1],xyzs[dof[2]-1],xyzs[dof[3]-1],xyzs[dof[4]-1]])
            mirror_vec[dof[1]-1]+=dof[5]*vectors[0]
            mirror_vec[dof[2]-1]+=dof[5]*vectors[1]
            mirror_vec[dof[3]-1]+=dof[5]*vectors[2]
            mirror_vec[dof[4]-1]+=dof[5]*vectors[3]

    for i in range(nAthoms):
        v_len=np.linalg.norm(mirror_vec[i])
        if v_len > 0.01:
            np.multiply(1/v_len, mirror_vec[i])


    mul_res=np.sum(mirror_vec*m_grad)
    sqr_res=np.sum(mirror_vec*mirror_vec)
    
    mirror_grad_cos=abs(mul_res/(sqr_res*np.sum(m_grad*m_grad))**0.5)
    mirror_grad_cos_used=mirror_grad_cos/5
    if(be_verbose):
        print(abs(mul_res/(sqr_res*np.sum(m_grad*m_grad))**0.5))
        print (f"mgcos {mirror_grad_cos_used}")
    m_grad=np.subtract(m_grad,(1+mirror_grad_cos_used)*np.multiply(mul_res/sqr_res,mirror_vec))#Это и будет эффективная отражнная сила
    m_grad_mean=np.array([0,0,0])#Вычтем "движение центра масс"
    for i in range(nAthoms):
        m_grad_mean=np.add(m_grad_mean,m_grad[i])
    m_grad_mean/=nAthoms
    
    for i in range(nAthoms):
        m_grad[i]-=m_grad_mean
    
    return m_grad, mirror_grad_cos

if __name__ == "__main__":
    '''#sn2
    xyzs=[
        [0.45620232128346, -0.99341949954807, -0.28183877557043],
        [0.85185551625923,  0.66887978618482,  0.85025682681636],
        [1.95347714750287, -2.70008895577006, -1.05030770331934]
    ]
    forces=[
        [ 2.5138235359462E-02,  -2.5999944410632E-02,  -7.1766747363902E-03],
        [ 6.7711430219807E-03,  -7.8534391528971E-03,  -4.5208342535619E-03],
        [-2.8717945882720E-02,   3.2351557728985E-02,   1.3758330914944E-02]        
    ]
    search_dofs=[
        ["b", 1, 2, 1],
        ["b", 2, 3, -1]
    ]
    print(mirror_fn(forces,xyzs,search_dofs))
    '''
    xyzs=[
        [-2.03398264799198,  0.09094558445260, -0.50127577437998],
        [ 0.95985948829336,  0.23016804578418, -0.60312180405917],
        [-1.69975909975752, -2.04599168483542,  0.75360343800011],
        [-0.42275103465649, -2.10827034835032,  1.08807203025491]
    ]
    forces=[
        [-1.2092892025177E-03,   7.7317952398507E-03,  -4.5450736506581E-03],
        [ 2.4532704085565E-04,   2.8422499997434E-04,  -2.3606650872549E-04],
        [ 1.2073287540207E-03,  -7.7143588466825E-03,   4.5259728780435E-03],
        [-1.9277947872847E-04,  -3.5061953151426E-04,   2.4079761792854E-04]
    ]
    search_dofs=[
        ["b", 1, 2, 1],
        ["b", 3, 4, 1]
    ]
    print(mirror_fn(forces,xyzs,search_dofs))
    