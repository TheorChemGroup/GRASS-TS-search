import numpy as np
forces=[
    [-2.5138235359462E-02,   2.5999944410632E-02,   7.1766747363902E-03],
    [-6.7711430219807E-03,   7.8534391528971E-03,   4.5208342535619E-03],
    [ 2.8717945882720E-02,  -3.2351557728985E-02,  -1.3758330914944E-02]
]
forces=-np.array(forces)
atoms=[
    [ 0.86209744464184,     -1.87729078100335,     -0.53259809712648],
    [-1.60977362347710,      1.26399960618399,      1.60675253827175],
    [ 3.69153679948519,     -5.10242863851781,     -1.98479390585265]
]
atoms=np.array(atoms)

dirCO=atoms[1]-atoms[0]
dirCCl=atoms[2]-atoms[0]
fCO=forces[1]-forces[0]
fCCl=forces[2]-forces[0]

faCO=np.dot(fCO,dirCO)/np.dot(dirCO,dirCO)
faCCl=np.dot(fCCl,dirCCl)/np.dot(dirCCl,dirCCl)

print(np.array(faCO, faCCl))
