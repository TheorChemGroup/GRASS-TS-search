import os, subprocess, time
works=["fullerene4_test"]

Ns=[]
gradtimes=[]
hesstimes=[]

cwd=os.getcwd()
for i in range(len(works)):
    path_folder=os.path.join(cwd,"tests",works[i])
    xyz_path=os.path.join(path_folder, "to_opt.xyz")
    with open(xyz_path,"r") as file:
        Ns.append(file.readline()[:-1])

    print(f"5hess {xyz_path}")
    time_start=time.time()
    for j in range(5):
        subprocess.call(["xtb", xyz_path, "--hess"],stdout=None)
    time_end=time.time()
    with open(os.path.join(path_folder,"hess_calc"),"w+") as file:
        file.write(str((time_end-time_start)/5))
    hesstimes.append(str((time_end-time_start)/5))
    
    print(f"100grad {xyz_path}")
    time_start=time.time()
    for j in range(100):
        subprocess.call(["xtb", xyz_path, "--grad"],stdout=None)
    time_end=time.time()
    with open(os.path.join(path_folder,"grad_calc"),"w+") as file:
        file.write(str((time_end-time_start)/100))
    gradtimes.append(str((time_end-time_start)/100))
    
with open("table_gradhess_f4.csv","w+") as file:
    file.write("N, gradtime, hesstime\n")
    for i in range(len(Ns)):
        file.write(f"{Ns[i]}, {gradtimes[i]}, {hesstimes[i]}\n")
