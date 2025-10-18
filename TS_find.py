
import numpy as np,os,subprocess,datetime,copy
from pathlib import Path
class usingMethod:
    def __init__(self,programm_name, dict_of_pars):
        self.programm_name=programm_name
        self.settings=dict()

        if self.programm_name=="xtb":
            self.settings["chrg"]=dict_of_pars["chrg"]
            self.settings["uhf"]=dict_of_pars["uhf"]
            self.settings["force_constant"]=dict_of_pars["force_constant"]
            self.settings["nAtoms"]=dict_of_pars["nAtoms"]
            self.settings["solvent"]=dict_of_pars["solvent"]
            self.settings["rpath"]=dict_of_pars["rpath"]
            self.settings["acc"]=dict_of_pars["acc"]
        elif self.programm_name=="orca":
            self.settings["method_str"]=dict_of_pars["method_str"]#i.e. "B3LYP def2-SVP"
            self.settings["memory"]=dict_of_pars["memory"]
            self.settings["nprocs"]=dict_of_pars["nprocs"]
            self.settings["chrg"]=dict_of_pars["chrg"]
            self.settings["mult"]=dict_of_pars["mult"]
            self.settings["nAtoms"]=dict_of_pars["nAtoms"]
            self.settings["solvent"]=dict_of_pars["solvent"]
            self.settings["rpath"]=dict_of_pars["rpath"]
            self.ORCA_PATH=dict_of_pars["ORCA_PATH"]
        else:
            raise ValueError(f"Unknown method: {self.programm_name}")
        
    #generic functions
    def get_energy(self):
        if self.programm_name=="xtb":
            return float(self.xyzs_strs[1].split()[1])
        if self.programm_name=="orca":
            return self.__get_energy_orca()
    def read_xyz(self, xyz_name):
        if self.programm_name=="xtb":
            if xyz_name=="!result":
                xyz_name="xtbopt.xyz"
            self.xyzs_strs=self.__read_file(xyz_name)
        if self.programm_name=="orca":
            if xyz_name=="!result":
                xyz_name="inpfile_trj.xyz"
                self.xyzs_strs=self.__read_last_struct(xyz_name)
            else:
                self.xyzs_strs=self.__read_file(xyz_name)

    def read_grad(self):
        if self.programm_name=="xtb":
            self.grad_strs=self.__read_file("gradient")
        elif self.programm_name=="orca":
            self.grad_strs=self.__read_file("inpfile.engrad")

    def extract_AB_dir(self, num_A:int,num_B:int):
        if self.programm_name=="xtb":
            return self.__extract_AB_dir_xtb(num_A,num_B)
        elif self.programm_name=="orca":
            return self.__extract_AB_dir_orca(num_A,num_B)
    
    def angle_3_ath(self, a1, a2, a3):
        v1=self.extract_AB_dir(a2,a1)
        v2=self.extract_AB_dir(a2,a3)
        cosval=np.matmul(v1,np.array(v2).T)/ (np.linalg.norm(v1)*np.linalg.norm(v2))
        if cosval>1:
            return 0
        elif cosval<-1:
            return np.pi
        else:
            return np.arccos(cosval)
    
    def d_4_ath(self, a1, a2, a3, a4): #from wikipedia
        u1=self.extract_AB_dir(a2,a1)
        u2=self.extract_AB_dir(a3,a2)
        u3=self.extract_AB_dir(a4,a3)
        u1pu2=np.cross(u1,u2)
        u2pu3=np.cross(u2,u3)
        arg1=np.matmul(u2,np.cross(u1pu2,u2pu3).T)
        arg2=np.linalg.norm(u2)*np.matmul(u1pu2,u2pu3.T)
        return -np.arctan2(arg1, arg2)
        
        
    def extractGradient(self, num):
        if self.programm_name=="xtb":
            return self.__extractGradient_xtb(num)
        elif self.programm_name=="orca":
            return self.__extractGradient_orca(num)
        
    def grad(self,xyz_name):
        if self.programm_name=="xtb":
            if xyz_name=="!result":
                xyz_name="xtbopt.xyz"
            self.__grad_xtb(xyz_name)
        elif self.programm_name=="orca":
            return #градиент выносится в файл в оптимизации в opt_constrain 
            self.__grad_orca(xyz_name)


    def opt_constrain(self,xyz_name:str,constrains:list):
        '''constrains is list of same lists: [ctype("bond","angle"),[related athoms],value]'''
        if self.programm_name=="xtb":
            if xyz_name=="!result":
                xyz_name="xtbopt.xyz"
            control_strs=[]
            for constrain in constrains:
                if constrain[0]=="bond":
                    control_strs.append(f"    distance: {constrain[1][0]}, {constrain[1][1]}, {constrain[2]}\n")
                elif constrain[0]=="angle":
                    control_strs.append(f"    angle: {constrain[1][0]}, {constrain[1][1]}, {constrain[1][2]}, {constrain[2]}\n")
                elif constrain[0]=="dihedral":
                    control_strs.append(f"    dihedral: {constrain[1][0]}, {constrain[1][1]}, {constrain[1][2]}, {constrain[1][3]}, {constrain[2]}\n")
                else:
                    raise ValueError(f"Unknown constrain type: {constrain[0]}")
            self.__save_control_xtb(control_strs)
            self.__opt_xtb(xyz_name)
        elif self.programm_name=="orca":
            self.__opt_constrain_orca(xyz_name,constrains)

    def xmol_xyzs_strs(self):
        if self.programm_name=="xtb":
            return self.xyzs_strs
        if self.programm_name=="orca":
            return self.xyzs_strs
        

    def __read_file(self,file_name:str):
        file_strs=[]
        with open(os.path.join(self.settings["rpath"],file_name),"r") as file:
            line=file.readline()
            while line!="":
                file_strs.append(line)
                line=file.readline()
        return file_strs
    
    def __read_last_struct(self,xyz_name):
        file_strs=self.__read_file(xyz_name)
        return file_strs[(-2-self.settings["nAtoms"]) : ]
    #~generic functions

    #xtb
    def __save_control_xtb(self,control_strs):
        with open(os.path.join(self.settings["rpath"],"control"),"w+") as control:
            control.writelines([f'$chrg {self.settings["chrg"]}\n',"$constrain\n"])
            control.writelines(control_strs)
            control.writelines([f'    force constant = {self.settings["force_constant"]}\n',"$end\n"])
    
    def __opt_xtb(self,xyz_name):
        with open(os.path.join(self.settings["rpath"],"xtbout"),"w+") as xtbout:
            if self.settings["solvent"]=="vacuum":
                p=subprocess.call(["xtb", "gfn1", xyz_name, "-I", "control","--uhf", str(self.settings["uhf"]),"--acc", str(self.settings["acc"]), "--opt"],stdout=xtbout)
            else:
                p=subprocess.call(["xtb", "gfn1", xyz_name, "-I", "control","--uhf", str(self.settings["uhf"]),"--alpb",self.settings["solvent"],"--acc", str(self.settings["acc"]), "--opt"],stdout=xtbout)
            if p!=0:
                print("abnormal termination of xtb. Exiting")
                os.chdir(self.initial_cwd)
                raise(Exception)
    def __grad_xtb(self,xyz_name):
        with open(os.path.join(self.settings["rpath"],"xtbout"),"w+") as xtbout:
            if self.settings["solvent"]=="vacuum":
                p=subprocess.call(["xtb", "gfn1", xyz_name, "--chrg", str(self.settings["chrg"]), "--uhf", str(self.settings["uhf"]),"--acc", str(self.settings["acc"]),"--grad"],stdout=xtbout)
            else:
                p=subprocess.call(["xtb", "gfn1", xyz_name, "--chrg", str(self.settings["chrg"]), "--uhf", str(self.settings["uhf"]),"--alpb", self.settings["solvent"],"--acc", str(self.settings["acc"]),"--grad"],stdout=xtbout)
            if p!=0:
                print("abnormal termination of xtb. Exiting")
                os.chdir(self.initial_cwd)
                raise(Exception)


    def __extractGradient_xtb(self, num):
        line_num=num+self.settings["nAtoms"]+1
        gradline=self.grad_strs[line_num]
        gradstr_arr=gradline.split()
        gradarr=[]
        for i in range(len(gradstr_arr)):
            gradarr.append(float(gradstr_arr[i]))
        return np.array(gradarr)

    def __extract_AB_dir_xtb(self, num_A:int,num_B:int):
        vec_A=self.xyzs_strs[num_A+1].split()[1:]
        for num,coord in enumerate(vec_A):
            vec_A[num]=float(coord)
    
        vec_B=self.xyzs_strs[num_B+1].split()[1:]
        for num,coord in enumerate(vec_B):
            vec_B[num]=float(coord)
    
        res=np.subtract(vec_B,vec_A)
        return res
    #~xtb
    #orca
    def __get_energy_orca(self):
        return float(self.grad_strs[7])#в 8-й строке находится значение энергии
    
    def __extract_AB_dir_orca(self, num_A, num_B):
        vec_A=self.xyzs_strs[num_A+1].split()[1:]
        for num,coord in enumerate(vec_A):
            vec_A[num]=float(coord)
    
        vec_B=self.xyzs_strs[num_B+1].split()[1:]
        for num,coord in enumerate(vec_B):
            vec_B[num]=float(coord)
    
        res=np.subtract(vec_B,vec_A)
        return res
    
    def __extractGradient_orca(self, num):
        grad_begin=11#12-я строка - X координата градиента первого атома
        grad_vec=[]
        for i in range(grad_begin+(num-1)*3, grad_begin+(num)*3):
            grad_vec.append(float(self.grad_strs[i]))
        return np.array(grad_vec)
    
    def __grad_orca(self,xyz_name):
        self.__make_and_run_orca_job(xyz_name,"grad")
    def __opt_constrain_orca(self,xyz_name,constrains):
        self.__make_and_run_orca_job(xyz_name,"opt", constrains=constrains)
        
    def __make_and_run_orca_job(self,xyz_name,target,constrains=[], jobname="inpfile.inp"):
        if target=="opt":
            target_str="Opt"
        elif target=="grad":
            target_str="EnGrad"
        elif target=="TS":
            target_str="OptTS"
        else:
            raise ValueError(f"unknown target: {target}")
        
        job=[f'%maxcore {self.settings["memory"]} %pal nprocs {self.settings["nprocs"]} end\n',
            f'! {self.settings["method_str"]} {"" if self.settings["solvent"]=="vacuum" else "CPCM("+self.settings["solvent"]+")"} {target_str}\n']
        if constrains!=[]:
            job.extend(["%geom\n", "Constraints\n"])
            for constrain in constrains:
                if constrain[0]=="bond":
                    job.append("{B "+f"{constrain[1][0]-1} {constrain[1][1]-1} {constrain[2]}" +" C}\n")
                elif constrain[0]=="angle":
                    job.append("{A "+f"{constrain[1][0]-1} {constrain[1][1]-1} {constrain[1][2]-1} {constrain[2]}" +" C}\n")
                elif constrain[0]=="dihedral":
                    job.append("{D "+f"{constrain[1][0]-1} {constrain[1][1]-1} {constrain[1][2]-1} {constrain[1][3]-1} {constrain[2]}" +" C}\n")
                else:
                    raise ValueError(f"Unknown constrain type: {constrain[0]}")
            job.extend(["end\n","end\n"])
        if xyz_name=="!result":
            xyzs_in=self.xyzs_strs[2:]
        else:
            xyzs_in=self.__read_file(xyz_name)
            if xyzs_in[0].split()[0].isdigit():#xmol
                xyzs_in=xyzs_in[2:]
        job.append(f'*xyz {self.settings["chrg"]} 1\n')
        job.extend(xyzs_in)
        job.extend(["\n","*\n"])
        with open(jobname, "w+") as file:
            file.writelines(job)
        with open(os.path.join(self.settings["rpath"],"outfile.out"),"w+") as orcaout:
            subprocess.call([os.path.join(self.ORCA_PATH,"orca"), jobname,"--use-hwthread-cpus --oversubscribe"],stdout=orcaout)    

    #~orca
    
class optTS:
    def __init__(self, xyz_path:str,threshold_force:float=0, threshold_rel:float=0,programm=dict(name="xtb"), maxstep:int=7000, print_output:bool=True, mode="strict"):
        cwd=os.getcwd()
        
        rpath=os.path.join(cwd, os.path.dirname(xyz_path))
        xyz_name=os.path.basename(xyz_path)

        
        if threshold_force==0 and threshold_rel==0:
            print("please, enter threshold_force or (and) threshold_rel")
            return
        self.const_settings=dict(rpath=rpath,xyz_name=xyz_name,print_output=print_output, threshold_force=threshold_force,threshold_rel=threshold_rel,maxstep=int(maxstep), mode=mode)
        self.settings=dict(step=0,prev_dc=100,bond_reach_critical_len=True)

        #self.log("",os.path.join(self.const_settings["rpath"],"log_doc"))
        #self.log("",os.path.join(self.const_settings["rpath"],"log_E"))
        #self.log("",os.path.join(self.const_settings["rpath"],"way_log.txt"))

        self.initial_cwd = os.getcwd()
        os.chdir(self.const_settings["rpath"])
    
        self.logname=os.path.join(self.initial_cwd,"grad_log")
    
        with open(self.logname,"w+") as file:
            file.write( (datetime.datetime.now()).strftime("%Y m%m d%d %H:%M:%S\n") )
    
        self.search_DoFs=[] # список СС, для которых ищется ПС
        self.ifprint("reading inputs\n")
    
        self.read_DoFs()
        self.const_settings["nDoFs"]=len(self.search_DoFs)
        if self.const_settings["nDoFs"]<3 and self.const_settings["mode"]=="autostrict":
            self.const_settings["mode"]="strict"

        self.log_xyz("new") 
        
        #Все длины, над которыми производятся операции - в ангстремах
        self.init_DoFs={}#[(a, b)] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
        self.lens={}#[(a, b)] текущие длины связей (к которым применяется изменение длины по градиенту)
    
        self.xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
        self.xyzs_strs=self.read_file(self.const_settings["xyz_name"])
        self.const_settings["nAtoms"]=len(self.xyzs_strs)-2
        
        if programm["name"]=="xtb":
            dict_to_uM=dict(rpath=self.const_settings["rpath"],
                            chrg=self.const_settings["chrg"],
                            uhf=self.const_settings["mult"]-1,
                            force_constant=programm["force_constant"],
                            acc=programm["acc"],
                            solvent=self.const_settings["solvent"],
                            nAtoms=self.const_settings["nAtoms"])
        elif programm["name"]=="orca":
            dict_to_uM=dict(method_str=programm["method_str"],
                            memory=programm["memory"],
                            nprocs=programm["nprocs"],
                            ORCA_PATH=programm["ORCA_PATH"],
                            
                            rpath=self.const_settings["rpath"],
                            chrg=self.const_settings["chrg"],
                            mult=self.const_settings["mult"],
                            solvent=self.const_settings["solvent"],
                            nAtoms=self.const_settings["nAtoms"])
        else:
            raise ValueError(f'Unknown method {programm["name"]}')

        self.Method=usingMethod(programm["name"],dict_to_uM)
        self.Method.read_xyz(self.const_settings["xyz_name"])
        self.log_xyz() 

        self.change_projections={}
        self.change_projections["vn_for_change"]=[np.array([self.const_settings["nDoFs"]**-0.5 for i in range(self.const_settings["nDoFs"])])]#self.make_ort111(self.const_settings["nDoFs"])
        for i in range(self.const_settings["nDoFs"]-1):
            self.change_projections["vn_for_change"].append([0 for j in range(self.const_settings["nDoFs"])])
        self.change_projections["reliability"]=[(2 if i>0 else 10) for i in range(self.const_settings["nDoFs"])]
        self.change_projections["selected_vector"]=1
        
        self.init_DoFs=self.find_reac_type_by_phases__and__measure_init_DoFs()
        self.init_change_projections=copy.deepcopy(self.change_projections)
        

        if self.const_settings["threshold_force"]=="auto":
            self.const_settings["threshold_force"]=self.mean_force()
            self.ifprint(f'because optimized cap is \"auto\", calculated threshold_force is {self.const_settings["threshold_force"]}')

        self.not_completed=True
        self.proceed()
    
    #mathematics
    @staticmethod  
    def projection(va,vb):
        return np.multiply( np.matmul(va,vb.T)/np.matmul(vb,vb.T), vb )#a to b projection
    @staticmethod  
    def vec_len(v):
        return (np.matmul(v,v.T))**0.5
    @staticmethod               
    def sign(x:float):
        if x>0:
            return 1
        if x<0:
            return -1
        return 0
    @staticmethod
    def vsign(v1:list,v2:list):
        cos_v1v2=np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
        if cos_v1v2>0.9:
            return 1
        if cos_v1v2<-0.9:
            return -1
        return None
        
    @staticmethod               
    def cos2v(v1:list, v2:list):
        return abs(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    
    def change_fn(self,length:float, cap:float):    
        return min((length/self.const_settings["nDoFs"]+10*length**2+1000*length**3)*2, cap)
        
    
    def produce_new_vector(self,changes):#процесс Грамма-Шмидта
        init=copy.deepcopy(changes)
        self.ifprint(init)
        for vector in self.change_projections["vn_for_change"]:
            if np.linalg.norm(vector)>0.1:
                init=np.subtract(init,self.projection(np.array(init),np.array(vector)))
            else:
                break
        init=np.multiply(1/np.linalg.norm(init),init)
        self.ifprint(init)
        return init
    #~mathematics
 
    #rw
    def ifprint(self,to_print):
        if self.const_settings["print_output"]:
            print(to_print)

    def log_xyz(self, mode=None):
        if mode=="new":
            with open(os.path.join(self.const_settings["rpath"],"TS_search_log.xyz"),"w+") as log:
                pass
        else:
            with open(os.path.join(self.const_settings["rpath"],"TS_search_log.xyz"),"a") as file:
                print(os.path.join(self.const_settings["rpath"],"TS_search_log.xyz"))
                file.writelines(self.Method.xmol_xyzs_strs())
    
    @staticmethod
    def log(str:str,logname:str):
        with open(logname,"a" if str!="" else "w+") as file:
            file.write(str)
    
    def read_file(self,file_name:str):
        file_strs=[]
        with open(os.path.join(self.const_settings["rpath"],file_name),"r") as file:
            line=file.readline()
            while line!="":
                file_strs.append(line)
                line=file.readline()
        return file_strs
    #~rw

    #init fns
    def read_DoFs(self):
        self.search_DoFs=[]
        with open(os.path.join(self.const_settings["rpath"],"bonds_to_search"),"r") as bonds:
            self.const_settings["chrg"]=int(bonds.readline())
            multline=bonds.readline()
            if(multline.startswith("auto")):
                self.const_settings["mult"]=abs(self.const_settings["chrg"])+1
            else:    
                self.const_settings["mult"]=int(multline)
            self.const_settings["solvent"]=bonds.readline().split()[0]
            line=bonds.readline()
            while line != "":
                if line.startswith("b") or line.startswith("a") or line.startswith("d"):#это СС
                    line_split=line.split()
                    if line_split[0]=='b':
                        if int(line_split[1])<int(line_split[2]):
                            self.search_DoFs.append(["b", int(line_split[1]), int(line_split[2]), float(line_split[3])])
                        else:
                            self.search_DoFs.append(["b", int(line_split[2]), int(line_split[1]), float(line_split[3])])
                    elif line_split[0]=='a':
                        if int(line_split[1])<int(line_split[3]):
                            self.search_DoFs.append(["a", int(line_split[1]), int(line_split[2]), int(line_split[3]), float(line_split[4])])
                        else:
                            self.search_DoFs.append(["a", int(line_split[3]), int(line_split[2]), int(line_split[1]), float(line_split[4])])
                    elif line_split[0]=='d':
                        if int(line_split[1])<int(line_split[4]):
                            self.search_DoFs.append(["d", int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]), float(line_split[5])])
                        else:
                            self.search_DoFs.append(["d", int(line_split[4]), int(line_split[3]), int(line_split[2]), int(line_split[1]), float(line_split[5])])
                
                line=bonds.readline()
    
    def find_reac_type_by_phases__and__measure_init_DoFs(self):#именно то, что написано на упаковке, более короткого, но осмысленного названия я придумать не смог
        RAD2DEG=57.295779513
        init_DoFs={}
        reac_type=2
        phases_vec=[]
        for DoF in self.search_DoFs:
            if DoF[0]=="b":
                phases_vec.append(DoF[3])
                key = (DoF[1], DoF[2])   
                init_DoFs[key]=self.vec_len(self.Method.extract_AB_dir(DoF[1],DoF[2])) 
            elif DoF[0]=="a":
                phases_vec.append(DoF[4])
                key = (DoF[1], DoF[2], DoF[3])   
                init_DoFs[key]=self.Method.angle_3_ath(DoF[1], DoF[2], DoF[3])*RAD2DEG
            elif DoF[0]=="d":
                phases_vec.append(DoF[5])
                key = (DoF[1], DoF[2], DoF[3], DoF[4])   
                init_DoFs[key]=self.Method.d_4_ath(DoF[1], DoF[2], DoF[3], DoF[4])*RAD2DEG

        print(init_DoFs)

        
        self.phases_vec=np.multiply(1/np.linalg.norm(phases_vec),phases_vec)
        self.ifprint(self.phases_vec)
        self.ifprint(init_DoFs)
        return init_DoFs
    #~init fns

    #reliability
    def increase_rel(self,vec_num:int):#увеличиваем уверенность в том, что это правильное изменение по знаку
        REL_MAX=60
        REL_INIT=2
        if self.change_projections["reliability"][vec_num]<=1:
            self.change_projections["reliability"][vec_num]=REL_INIT
        else:
            self.change_projections["reliability"][vec_num]+=1
            if self.change_projections["reliability"][vec_num]>REL_MAX:
                self.change_projections["reliability"][vec_num]=REL_MAX
    def decrease_rel(self,vec_num:int):#уменьшаем уверенность в том, что это правильное изменение по знаку
        if self.change_projections["reliability"][vec_num]<=1:
            if self.const_settings["mode"]!="strict":
                self.change_projections["signs"][vec_num]=-self.change_projections["signs"][vec_num]
                self.increase_rel(vec_num)
        else:
            if self.change_projections["reliability"][vec_num] >= 12:
                self.change_projections["reliability"][vec_num]-=2
            else:    
                self.change_projections["reliability"][vec_num]-=1
            if self.const_settings["nDoFs"]==1:
                self.change_projections["reliability"][vec_num]-=1
    #~reliability
        
    #main loop fns
    def reset(self):
        self.settings["prev_dc"]=100
        self.settings["bond_reach_critical_len"]=True
        self.change_projections=copy.deepcopy(self.init_change_projections)

                
    def proceed(self):
        while self.not_completed:
            if self.settings["bond_reach_critical_len"]==True:
                self.lens.clear()
                self.ifprint("lens is clear")
    
                self.constrain_list=[]
                for DoF_atoms, DoF_value in zip(self.init_DoFs.keys(), self.init_DoFs.values()):
                    if len(DoF_atoms)==2:
                        const_type="bond"
                    elif len(DoF_atoms)==3:
                        const_type="angle"
                    elif len(DoF_atoms)==4:
                        const_type="dihedral"               
                    self.constrain_list.append([const_type,DoF_atoms, DoF_value])

                self.Method.opt_constrain(self.const_settings["xyz_name"],self.constrain_list)
                self.Method.read_xyz("!result")
                self.log_xyz()
                
                self.settings["bond_reach_critical_len"]=False
    
            else:
                

                #self.log(f"{self.xyzs_strs[1].split()[1]}\n",os.path.join(self.const_settings["rpath"],"log_E"))
        
                if self.const_settings["nDoFs"]>1:
                    string_curve=f'{self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[0][1],self.search_DoFs[0][2]))} {self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[1][1],self.search_DoFs[1][2]))} {self.Method.get_energy()}\r\n'
                    self.log(string_curve, "way_log.txt")
                
                self.constrain_list=[]
                
                proj_len=self.move_DoFs()
                self.not_completed = not self.check_tresholds_converged(proj_len)
                
                if self.settings["step"]>=self.const_settings["maxstep"]:
                    self.ifprint("\033[91mnot optimized but reached maximum number of steps\033[00m") 
                    self.not_completed=False
                
                if self.not_completed:
                    self.ifprint("opt geometry with new constrains") 
                    self.Method.opt_constrain("!result",self.constrain_list)
                    self.Method.read_xyz("!result")
                    self.log_xyz()

            if self.not_completed:
                self.ifprint(f'\nstep {self.settings["step"]}')
                self.ifprint("gradient calculation")
            
                self.Method.grad("!result")
                self.Method.read_grad() 
            else:
                self.log(f"completed at {(datetime.datetime.now()).strftime('%Y m%m d%d %H:%M:%S')}\n",self.logname)
    
            
        os.chdir(self.initial_cwd)

    def move_DoFs(self):
        MIN_BOND=0.8
        MAX_BOND=3.5
        KOEF_ELEM_ANGLE_ACT=8#коэффициент перевода из углового действия [рад] в элементарное [Ангст]
        KOEF_ELEM_DIHEDRAL_ACT=15
        MIN_ANG=0.4#[rad]
        RAD2DEG=57.295779513
        sum_changes=0
        min_change=1000
        max_change=-1000
        changes=[]
        self.settings["step"]+=1
        for i,DoF in enumerate(self.search_DoFs):#найдём градиент (желание растянуться) вдоль каждой связи
            if DoF[0]=="b":
                num_A=DoF[1]
                num_B=DoF[2]
                key=(num_A, num_B)
    
                grad_A=self.Method.extractGradient(num_A)
                grad_B=self.Method.extractGradient(num_B)
                AB_dir=self.Method.extract_AB_dir(num_A,num_B)
                summ_grad=np.subtract(grad_A,grad_B)
    
                if not key in self.lens.keys():
                    self.ifprint(f"key \"{key}\" not in lens.keys()") 
                    self.lens[key] = self.vec_len(AB_dir)
    
                s_g_proj=self.projection(summ_grad, AB_dir)
                proj_len=self.vec_len(s_g_proj)
    
                s_g_p_sign=self.vsign(s_g_proj,AB_dir)
    
                changes.append(s_g_p_sign*proj_len)#удлиннение (если отрицательно - укорочение) связи
                max_change=max(changes[i], max_change)
                min_change=min(changes[i], min_change)
                sum_changes+=changes[i]
            elif DoF[0]=="a":#Здесь необходимо привести всё к некоторому "элементарному дейтвию" - величине, отражающей смещение вдоль каждой координаты в равной степени, как растяжения, так и повороты
                num_A, num_B, num_C = DoF[1],DoF[2],DoF[3]
                key = (num_A, num_B, num_C)
                if not key in self.lens.keys():
                    self.ifprint(f"key \"{key}\" not in lens.keys()") 
                    self.lens[key] = self.Method.angle_3_ath(num_A, num_B, num_C)*RAD2DEG

                grad_A=self.Method.extractGradient(num_A)
                grad_C=self.Method.extractGradient(num_C)
                BA_dir=self.Method.extract_AB_dir(num_B,num_A)
                BC_dir=self.Method.extract_AB_dir(num_B,num_C)
                BA_len=np.linalg.norm(BA_dir)
                BC_len=np.linalg.norm(BC_dir)
                norm_ABC=np.cross(BA_dir,BC_dir)
                norm_ABC=np.multiply(1/np.linalg.norm(norm_ABC), norm_ABC)

                force_A=np.subtract(self.projection(grad_A, norm_ABC), grad_A)
                force_C=np.subtract(self.projection(grad_C, norm_ABC), grad_C)
                m_F_A=np.cross(BA_dir,force_A)
                m_F_C=np.cross(BC_dir,force_C)
                sum_m=np.subtract(np.multiply(1/BC_len,m_F_C), np.multiply(1/BA_len,m_F_A))
                sum_m_len=np.linalg.norm(sum_m)
                s_m_sign=self.vsign(sum_m, norm_ABC)
                changes.append(s_m_sign*sum_m_len*KOEF_ELEM_ANGLE_ACT)
                max_change=max(changes[i], max_change)
                min_change=min(changes[i], min_change)
            elif DoF[0]=="d":#Здесь необходимо привести всё к некоторому "элементарному дейтвию" - величине, отражающей смещение вдоль каждой координаты в равной степени, как растяжения, так и повороты
                num_A, num_B, num_C, num_D = DoF[1],DoF[2],DoF[3],DoF[4]
                key = (num_A, num_B, num_C, num_D)
                if not key in self.lens.keys():
                    self.ifprint(f"key \"{key}\" not in lens.keys()") 
                    self.lens[key] = self.Method.d_4_ath(num_A, num_B, num_C, num_D)*RAD2DEG

                grad_A=self.Method.extractGradient(num_A)
                grad_D=self.Method.extractGradient(num_D)
                BA_dir=self.Method.extract_AB_dir(num_B,num_A)
                BC_dir=self.Method.extract_AB_dir(num_B,num_C)
                CD_dir=self.Method.extract_AB_dir(num_C,num_D)

                BAr_dir=np.subtract(BA_dir, self.projection(BA_dir, BC_dir))
                CDr_dir=np.subtract(CD_dir, self.projection(CD_dir, BC_dir))
                BAr_len=np.linalg.norm(BAr_dir)
                CDr_len=np.linalg.norm(CDr_dir)
                
                force_A=np.subtract(self.projection(grad_A, BC_dir), grad_A)
                force_D=np.subtract(self.projection(grad_D, BC_dir), grad_D)
                m_F_A=np.cross(BAr_dir,force_A)
                m_F_D=np.cross(CDr_dir,force_D)
                sum_m=np.subtract(np.multiply(1/CDr_len,m_F_D), np.multiply(1/BAr_len,m_F_A))
                sum_m_len=np.linalg.norm(sum_m)
                s_m_sign=self.vsign(sum_m, BC_dir)
                changes.append(s_m_sign*sum_m_len*KOEF_ELEM_DIHEDRAL_ACT)
                max_change=max(changes[i], max_change)
                min_change=min(changes[i], min_change)

        div_of_changes=max_change-min_change#отклонение 
    
        self.ifprint(f'div of forces: \033[01mcur\033[00m {"{:.8f}".format(div_of_changes)}    \033[01mprev\033[00m {"{:.8f}".format(self.settings["prev_dc"])} (\033{"[92mless" if div_of_changes<self.settings["prev_dc"] else "[091mhigher"}\033[00m)') 
        
        inv_chang_to_TS=np.subtract(changes,np.multiply(2,self.projection(changes,self.phases_vec)))
        len_ic=np.linalg.norm(inv_chang_to_TS)
        inv_chang_to_TS=np.multiply(self.change_fn(len_ic,0.03)/len_ic,inv_chang_to_TS)
        for i,DoF in enumerate(self.search_DoFs):
            if DoF[0]=="b":
                key=(DoF[1], DoF[2])
                bond_len=self.lens[key]
                #блок того, что надо сделать для каждой связи
                bond_change=inv_chang_to_TS[i]
            
                #~блок того, что надо сделать для каждой связи
            
                res_bond=bond_len+bond_change
                self.ifprint(f'\033[93m{key[0]}, {key[1]}\033[00m\tchange {"{:14.10f}".format(bond_change)}, res {"{:14.10f}".format(res_bond)}')
                self.constrain_list.append(["bond",key, res_bond])
                self.lens[key]=res_bond
    
                if res_bond>MAX_BOND or res_bond<MIN_BOND:
                    self.reset()
                    #если в результате поиска ПС связь порвалась или замкнулась, то  результат сбрасывается, а в начальной геометрии соответствующая связь немного (на 0,02) растягивается или укорачивается, чтобы притяжение или отталкивание было меньше
                    key=(DoF[1], DoF[2])
    
                    if res_bond>MAX_BOND:
                        self.init_DoFs[key]-=0.1
                    else:
                        self.init_DoFs[key]+=0.1
                    self.ifprint(self.init_DoFs)
            elif DoF[0]=="a":
                key=(DoF[1], DoF[2], DoF[3])
                angle_len=self.lens[key]
                angle_change=inv_chang_to_TS[i]*RAD2DEG
                
                res_angle=angle_len+angle_change
                self.ifprint(f'\033[93m{key[0]}, {key[1]}, {key[2]}\033[00m\tchange {"{:14.10f}".format(angle_change)}, res {"{:14.10f}".format(res_angle)}')
                self.constrain_list.append(["angle",key, res_angle])
                self.lens[key]=res_angle
            elif DoF[0]=="d":
                key=(DoF[1], DoF[2], DoF[3], DoF[4])
                angle_len=self.lens[key]
                angle_change=inv_chang_to_TS[i]*RAD2DEG
                
                res_angle=angle_len+angle_change
                self.ifprint(f'\033[93m{key[0]}, {key[1]}, {key[2]}, {key[3]}\033[00m\tchange {"{:14.10f}".format(angle_change)}, res {"{:14.10f}".format(res_angle)}')
                self.constrain_list.append(["dihedral",key, res_angle])
                self.lens[key]=res_angle

        self.settings["prev_dc"]=div_of_changes
            
        self.log(f"{div_of_changes}\n",os.path.join(self.const_settings["rpath"],"log_doc"))
        return max(abs(min_change),max_change)
    
    
        
    def mean_force(self):
        search_atoms=set()
        for DoF in self.search_DoFs:
            if DoF[0]=="b":
                search_atoms.add(DoF[1])
                search_atoms.add(DoF[2])
            elif DoF[0]=="a":
                search_atoms.add(DoF[1])
                search_atoms.add(DoF[3])
            elif DoF[0]=="d":
                search_atoms.add(DoF[1])
                search_atoms.add(DoF[4])
        sum_forces=0
        num_forces=0
        for i in range(1,self.const_settings["nAtoms"]+1):
            if i not in search_atoms:
                vec_force=self.Method.extractGradient(i)
                sum_forces+=self.vec_len(vec_force)
                num_forces+=1
        return sum_forces/num_forces if num_forces else 10000
    
    def check_tresholds_converged(self,proj_len:float):
        trashold_template=lambda name,cur,target,conv:f'{name} trashold {"{:15.7f}".format(cur)} of {"{:15.7f}".format(target)}: \033{"[92m" if conv else "[91mnot "}converged\033[00m'
        
        converged=True
        if self.const_settings["threshold_force"]!=0:
            cond=proj_len<=self.const_settings["threshold_force"]
            self.ifprint(trashold_template("force   ",proj_len,self.const_settings["threshold_force"],cond))
            converged &= cond
                
        if self.const_settings["threshold_rel"]!=0:
            mean_not_opt=self.mean_force()
            cur_threshold_rel=proj_len/mean_not_opt
            cond=cur_threshold_rel<=self.const_settings["threshold_rel"]
            self.ifprint(trashold_template("relative",cur_threshold_rel,self.const_settings["threshold_rel"],cond))
            converged &= cond
               
        return converged
    #~main loop fns
#------run------#
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description='Method for finding TS by targeted bonds. You only need store bonds_to_search and <name>.xyz files to directory/ and then call that programm', epilog="When using ORCA, it's need to export its folder to PATH, LD_LIBRARY_PATH. If using multiprocessoring (openmpi) it's need to export its folders lib/ to LD_LIBRARY_PATH and bin/ to PATH")
    parser.add_argument("xyz_path", type=str, help="xmol .xyz file with structure. File can be in any directory")
    parser.add_argument("-tf", "--treshold-force", type=float, default=0.00004,dest="threshold_force", help="that trashold is converged when max force on optimizing bonds less than its value. Default: 0.00004")
    parser.add_argument("-tr", "--treshold-rel", type=float, default=8.,dest="threshold_rel", help="that trashold is converged when max force on optimizing bonds divided by mean force on unconstrained bonds less then its value. Default: 8")
    parser.add_argument("-mc", "--mirror-coef", type=float, default=1.,dest="mirror_coef", help="The projection of the force at reflection of the longitudinal component relative to the phase vector is multiplied by this value. A decrease leads to a decrease in velocity, while an increase can cause oscillations near the transition state. Default: 1")
    parser.add_argument("--verbose",const=True, default=False,action='store_const', help="print output")
    parser.add_argument("-s", "--steps", type=int, default=2000, help="maximum number of steps that allowed to optimize TS. Default: 2000")
    parser.add_argument("-p", "--programm", default="xtb",help="programm that used for gradient calculation and constraint optimization. \"xtb\" or \"orca\". Default: \"xtb\"")
    parser.add_argument("-xfc","--xtb-force-consant",type=float,default=6.,dest="xfc",help="if using xtb that force constant is used in control file. Default: 6")
    parser.add_argument("-acc","--acc",type=float,default=0.05, dest="acc",help="if using xtb that acc is used. Default: 0.05")
    parser.add_argument("-oms","--orca-method-string", type=str, default="B3LYP def2-SVP",dest="method_str", help="method string on the top of orca file. Default: \"B3LYP def2-SVP\"")
    parser.add_argument("-OPATH", "--ORCA-PATH", type=str, default="", dest="OPATH",help="PATH of ORCA. Is necessary for multiprocessor calculations")
    parser.add_argument("-onp","--orca-number-processors", type=int, default=1,dest="nprocs", help="number of processors that using in ORCA work. Default: 1")
    parser.add_argument("-omm","--orca-memory",type=int, default=2000, dest="mem", help="memory per processor amount using in ORCA work. Default: 2000")
    
    args=parser.parse_args()
    
    if args.nprocs > 1 and args.OPATH=="":
        print("It's NECESSARY to set OPATH when using multithreading with ORCA")
        exit(1)
    optTS(
          args.xyz_path, 
          threshold_rel=args.threshold_rel, 
          threshold_force=args.threshold_force, 
          print_output=args.verbose,
          maxstep=args.steps, 
          programm=dict(name=args.programm, 
                        method_str=args.method_str,
                        force_constant=args.xfc,
                        acc=args.acc,
                        nprocs=args.nprocs, 
                        memory=args.mem, 
                        ORCA_PATH=args.OPATH))
    '''
    initial_cwd=os.getcwd()
    optTS(xyz_path=os.path.join("scan_opt_Sn2Cl_many_dirs","work7", "to_opt.xyz"), threshold_rel=8, threshold_force=0.00004, print_output=True,mode="strict", maxstep=10**3, programm=dict(name="xtb", force_constant= 6))
    
    '''