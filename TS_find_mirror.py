
import numpy as np,os,subprocess,datetime,copy
from mirror_fn import mirror_fn
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
            return float(self.grad_strs[1].split()[6])
        if self.programm_name=="orca":
            return self.__get_energy_orca()
    def read_xyz(self, xyz_name):
        if self.programm_name=="xtb":
            if xyz_name=="!result":
                xyz_name="xtbopt.xyz"
            self.xyzs_strs=self.__read_file(xyz_name, self.settings["nAtoms"]+2)
        if self.programm_name=="orca":
            if xyz_name=="!result":
                xyz_name="inpfile_trj.xyz"
                self.xyzs_strs=self.__read_last_struct(xyz_name)
            else:
                self.xyzs_strs=self.__read_file(xyz_name, self.settings["nAtoms"]+2)

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
        return np.arctan2(arg1, arg2)
        
        
    def extractGradient(self, num):
        if self.programm_name=="xtb":
            return self.__extractGradient_xtb(num)
        elif self.programm_name=="orca":
            return self.__extractGradient_orca(num)
        
    def grad(self,xyz_name):
        if self.programm_name=="xtb":
            if xyz_name=="!result":
                xyz_name="xtbopt.xyz"
                with open(os.path.join(self.settings["rpath"],"xtbopt.xyz"), "w+") as file:
                    file.writelines(self.xyzs_strs)
            self.__grad_xtb(xyz_name)

        elif self.programm_name=="orca":
            self.__grad_orca(xyz_name)


    def opt_constrain(self,xyz_name:str,constrains:list):
        '''constrains is list of same lists: [ctype("bond","angle"),[related atoms],value]'''
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
        

    def __read_file(self,file_name:str, strs_cap=None):
        file_strs=[]
        
        with open(os.path.join(self.settings["rpath"],file_name),"r") as file:
            line=file.readline()
            if type(strs_cap)==int:
                str_num=0
                while line!="" and str_num<strs_cap:
                    file_strs.append(line)
                    line=file.readline()
                    str_num+=1
            else:
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
                raise(Exception)
    def __grad_xtb(self,xyz_name):
        with open(os.path.join(self.settings["rpath"],"xtbout"),"w+") as xtbout:
            if self.settings["solvent"]=="vacuum":
                p=subprocess.call(["xtb", "gfn1", xyz_name, "--chrg", str(self.settings["chrg"]), "--uhf", str(self.settings["uhf"]),"--acc", str(self.settings["acc"]),"--grad"],stdout=xtbout)
            else:
                p=subprocess.call(["xtb", "gfn1", xyz_name, "--chrg", str(self.settings["chrg"]), "--uhf", str(self.settings["uhf"]),"--alpb", self.settings["solvent"],"--acc", str(self.settings["acc"]),"--grad"],stdout=xtbout)
            if p!=0:
                print("abnormal termination of xtb. Exiting")
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
        job.append(f'*xyz {self.settings["chrg"]} {self.settings["mult"]}\n')
        job.extend(xyzs_in)
        job.extend(["\n","*\n"])
        with open(jobname, "w+") as file:
            file.writelines(job)
        with open(os.path.join(self.settings["rpath"],"outfile.out"),"w+") as orcaout:
            subprocess.call([os.path.join(self.ORCA_PATH,"orca"), jobname,"--use-hwthread-cpus"],stdout=orcaout)    

    #~orca

        
class optTS:
    def __init__(self, xyz_path:str,threshold_force:float=0, threshold_rel:float=0, mirror_coef:float=1, programm=dict(name="xtb"), mult=1, maxstep:int=7000, do_preopt=True,step_along=0, print_output:bool=True):
        cwd=os.getcwd()
        
        rpath=os.path.join(cwd, os.path.dirname(xyz_path))
        xyz_name=os.path.basename(xyz_path)

        
        if threshold_force==0 and threshold_rel==0:
            print("please, enter threshold_force or (and) threshold_rel")
            return
        self.const_settings=dict(rpath=rpath,xyz_name=xyz_name,print_output=print_output, threshold_force=threshold_force,threshold_rel=threshold_rel,maxstep=int(maxstep),mult=mult, do_preopt=do_preopt,step_along=step_along)
        self.settings=dict(step=0,prev_dc=100,bond_reach_critical_len=True, mirror_coef=mirror_coef)

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

        self.log_xyz("new") 
        
        #Все длины, над которыми производятся операции - в ангстремах
        self.init_DoFs={}#[(a, b)] (при этом a<b - номера атомов) - начальные длины связей. При разваливании ПС изменяются так, чтобы меньше разваливались (растягиваются при критичном сжатии, сжимаются при критичном растяжении..)
        self.lens={}#[(a, b)] текущие длины связей (к которым применяется изменение длины по градиенту)
    
        self.xyzs_strs=[]#строки координат атомов (вида "A X Y Z\n", A - символ элемента, X,Y,Z - его координаты)
        self.xyzs_strs=self.read_file(self.const_settings["xyz_name"])
        self.const_settings["nAtoms"]=int(self.xyzs_strs[0])
        
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
        os.chdir(self.initial_cwd)
    
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
            with open(os.path.join(self.const_settings["rpath"],"TS_search_m_log.xyz"),"w+") as log:
                pass
        else:
            with open(os.path.join(self.const_settings["rpath"],"TS_search_m_log.xyz"),"a") as file:
                self.ifprint(os.path.join(self.const_settings["rpath"],"TS_search_m_log.xyz"))
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


        for i in range(len(phases_vec)-1):
            if phases_vec[i]*phases_vec[i+1]<0:
                reac_type=1
                break
        self.ifprint(reac_type)

        self.change_projections["signs"]=[]
        self.phases_vec=np.multiply(1/np.linalg.norm(phases_vec),phases_vec)
        '''if reac_type==2:#как Дильс-Альдер   
            self.change_projections["signs"].append(1)
            for i in range(1,len(phases_vec)):
                self.change_projections["signs"].append(-1)
            self.phases_vec=self.change_projections["vn_for_change"][0]
        elif reac_type==1:#как sn2 
            self.change_projections["vn_for_change"][1]=self.produce_new_vector(phases_vec)
            self.change_projections["signs"].append(-1)
            self.change_projections["signs"].append(1)
            for i in range(2,len(phases_vec)):
                self.change_projections["signs"].append(-1)
            self.phases_vec=np.multiply(1/np.linalg.norm(phases_vec),phases_vec)'''
        self.ifprint(self.phases_vec)
        self.ifprint(init_DoFs)
        return init_DoFs
    
    def move_along(self,DoF_value,DoF_type,i):
        print((DoF_value,i,self.phases_vec))
        DoF_value*=1+self.const_settings["step_along"]*self.phases_vec[i]
        print(DoF_value)
        if DoF_type=="angle":
            DoF_value=min(DoF_value,179.9)
        elif DoF_type=="dihedral":
            DoF_value=min(DoF_value,179.9)
            DoF_value=max(DoF_value,-179.9)
        return DoF_value
    #~init fns
        
    #main loop fns
    def reset(self):
        self.settings["prev_dc"]=100
        self.settings["bond_reach_critical_len"]=True
        self.change_projections=copy.deepcopy(self.init_change_projections)

                
    def proceed(self):
        while self.not_completed:
            if self.settings["bond_reach_critical_len"]==True:
                
                self.Optimizer=0
                self.lens.clear()
                self.ifprint("lens is clear")
                self.atoms,self.xyzs=self.get_xyzs()

                self.prev_maxgrad=100000
                self.coef_grad=0.7
                
                #for ADAM
                self.mt=np.zeros((self.const_settings["nAtoms"],3))
                self.vt=0
                self.mt_r=np.zeros((self.const_settings["nAtoms"],3))
                self.vt_r=0
                self.diff_mt=np.zeros((self.const_settings["nAtoms"],3))
                self.diff_vt=0

                self.least_force=10e100

                if self.const_settings["do_preopt"]:
                    self.constrain_list=[]
                    for DoF_atoms, DoF_value,i in zip(self.init_DoFs.keys(), self.init_DoFs.values(), range(len(self.init_DoFs))):
                        if len(DoF_atoms)==2:
                            const_type="bond"
                        elif len(DoF_atoms)==3:
                            const_type="angle"
                        elif len(DoF_atoms)==4:
                            const_type="dihedral"                    

                        self.constrain_list.append([const_type,DoF_atoms, self.move_along(DoF_value,const_type,i)])
                    self.Method.opt_constrain(self.const_settings["xyz_name"],self.constrain_list)
                    self.Method.read_xyz("!result")
                
                self.atoms, self.xyzs=self.get_xyzs()
                self.log_xyz()
                self.Method.grad("!result")
                self.Method.read_grad() 
                
                #self.log("","way_log.txt")
                #str_way=""
                #for DoF_atoms in self.init_DoFs.keys():
                #    if len(DoF_atoms)==2:
                 #       str_way+=f"{np.linalg.norm(self.Method.extract_AB_dir(DoF_atoms[0],DoF_atoms[1]))} "
                #str_way+=f"{self.Method.angle_3_ath(5,1,8)} "
                #str_way+="\n"
                #self.log(str_way,"way_log.txt")
                

                self.settings["bond_reach_critical_len"]=False
                if self.const_settings["nDoFs"]>1:
                    string_curve=f'{self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[0][1],self.search_DoFs[0][2]))} {self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[1][1],self.search_DoFs[1][2]))} {self.Method.get_energy()}\r\n'
                    self.log(string_curve, "way_log.txt")
    
            else:
            
                self.constrain_list=[]
                
                proj_len=self.move_DoFs()

                if self.const_settings["nDoFs"]>1:
                    string_curve=f'{self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[0][1],self.search_DoFs[0][2]))} {self.vec_len(self.Method.extract_AB_dir(self.search_DoFs[1][1],self.search_DoFs[1][2]))} {self.Method.get_energy()}\r\n'
                    self.log(string_curve, "way_log.txt")

                if proj_len<self.least_force:
                    print(self.least_force)
                    self.least_force=proj_len
                    self.best_xyzs_strs=copy.deepcopy(self.Method.xyzs_strs)
                self.not_completed = not self.check_tresholds_converged(proj_len)
                
                #str_way=""
                #for DoF_atoms in self.init_DoFs.keys():
                #    if len(DoF_atoms)==2:
                #        str_way+=f"{np.linalg.norm(self.Method.extract_AB_dir(DoF_atoms[0],DoF_atoms[1]))} "
                #str_way+=f"{self.Method.angle_3_ath(5,1,8)} "
                #str_way+="\n"
                #self.log(str_way,"way_log.txt")
                
                if self.settings["step"]>=self.const_settings["maxstep"]:
                    self.ifprint("\033[91mnot optimized but reached maximum number of steps. Least-gradient geometry saved\033[00m")
                    with open(os.path.join(self.const_settings["rpath"],"result.xyz"),"w+") as file:
                        file.writelines(self.best_xyzs_strs)
                    self.not_completed=False
                
                if self.not_completed:
                    self.log_xyz()

            if self.not_completed:
                self.ifprint(f'\nstep {self.settings["step"]}')

            else:
                self.log(f"completed at {(datetime.datetime.now()).strftime('%Y m%m d%d %H:%M:%S')}\n",self.logname)
                with open(os.path.join(self.const_settings["rpath"],"result.xyz"),"w+") as file:
                    file.writelines(self.best_xyzs_strs)
            
        os.chdir(self.initial_cwd)

    def get_xyzs(self):
        atoms=[]
        xyzs=[]
        for xyz_str in self.Method.xyzs_strs[2:]:
            linesplit=xyz_str.split()
            atoms.append(linesplit[0])
            xyzs.append([float(linesplit[1]), float(linesplit[2]), float(linesplit[3])])
        return atoms, xyzs
    
    def get_grad(self):
        grad=[]
        maxgrad=0
        for i in range(self.const_settings["nAtoms"]):
            grad.append(self.Method.extractGradient(i+1))
            maxgrad=max(maxgrad,np.linalg.norm(grad[i]))
        return grad, maxgrad
    
    def alter_grad(self):
        maxgrad=0
        mingrad=1e200
        
        for i in range(self.const_settings["nAtoms"]):
            g_norm=np.linalg.norm(self.grad[i])
            maxgrad=max(maxgrad,g_norm)
            mingrad=min(mingrad,g_norm)

        strange_constant=1-(mingrad/maxgrad)**0.1#Эта величина 0..1. Чтобы вначале, когда mingrad/maxgrad большой - все атомы шевелились примерно одинаково быстро, а потом, когда он уменьшается - с той скоростью, с которой должны
        strange_constant=strange_constant**1.3#В целом, она улучшает сходимость реакций между 2 большими молекулами, помогая им повернуться/занять "молекулярно(т.е. в масштабе больших фрагментов)-правильное" положение
        #Физический смыcл как таковой в целом отсутствует - только алгоритмический // strange constant is strongly related to strange magick
        self.ifprint(f"strangeC {strange_constant}")
        for i in range(self.const_settings["nAtoms"]):
            g_norm=np.linalg.norm(self.grad[i])
            self.grad[i]=np.multiply((maxgrad/g_norm)**(strange_constant),self.grad[i])#Вот здесь при большом (около 1) strange_constant все едут примерно со скоростью max_grad (важно по сути только общее направление движения частей системы), а когда структура уже близка к стационарной точке (Strange_constant~=0) атомы двигаются так, как должны по градиенту
        

    def apply_grad(self):
        
        for i in range(self.const_settings["nAtoms"]):
            for j in range(3):
                self.xyzs[i][j]-=self.grad[i][j]*self.coef_grad
    
    def update_xyzs_strs(self):
        self.Method.xyzs_strs=self.Method.xyzs_strs[:1]
        self.Method.xyzs_strs.append(str(self.Method.get_energy())+"\n")
        for i in range(self.const_settings["nAtoms"]):
            self.Method.xyzs_strs.append(f"{self.atoms[i]} {float(self.xyzs[i][0])} {float(self.xyzs[i][1])} {float(self.xyzs[i][2])}\n")

    def mirror(self):
        self.grad, mirror_grad_cos=mirror_fn(self.grad,self.xyzs,self.search_DoFs, self.const_settings["print_output"])
        return mirror_grad_cos
    def move_DoFs(self):
        b1=0.2
        b1_diff=0.7
        b2=0.995
        eps=1e-7
        eta=1.5e-2
        eta_diff=1e-2
        eta_r=3e-2

        self.grad, maxgrad=self.get_grad()
        mgcos=self.mirror()
        mgsin_sqr=(1-mgcos*mgcos)**0.5
        
        TRUST_RAD=0.1
        #if(maxgrad*self.coef_grad>TRUST_RAD):
        #    self.coef_grad=TRUST_RAD/maxgrad
            
        self.alter_grad()
        if 1:
            #no_ADAM (rotational correction)
            self.mt = b1*self.mt + (1-b1)*self.grad#*self.coef_grad
            self.vt = b2*self.vt + (1-b2)*np.sum(self.grad*self.grad)#*self.coef_grad**2
            mt_bias=1/(1-b1**(self.settings["step"]+1))*self.mt
            vt_bias=(1-b2**(self.settings["step"]+1))*self.vt
            
            if(self.settings["step"]>0):
                cur_diff_mt = -self.prev_grad+self.grad
                self.diff_mt = b1_diff*self.diff_mt + (1-b1_diff)*cur_diff_mt
                self.diff_vt = b2*self.diff_vt + (1-b2)*np.sum(cur_diff_mt*cur_diff_mt)
                diff_mt_bias=1/(1-b1_diff**(self.settings["step"]+1))*self.diff_mt
                diff_vt_bias=(1-b2**(self.settings["step"]+1))*self.diff_vt

                vec_chang=-eta*(vt_bias+eps)**(-0.5) * mt_bias - eta_diff * mgsin_sqr * (diff_vt_bias+eps)**(-0.5)*diff_mt_bias
            else:
                vec_chang=-eta*(vt_bias+eps)**(-0.5) * mt_bias

            norm_chang=np.linalg.norm(vec_chang)
            if(norm_chang>TRUST_RAD):
                vec_chang=TRUST_RAD/norm_chang*vec_chang
            self.xyzs=self.xyzs+vec_chang
            
            self.prev_grad=copy.deepcopy(self.grad)
            
            self.update_xyzs_strs()
            self.Method.grad("!result")
            self.Method.read_grad() 
        elif 0:
            #ADAM (Egor idea)
            self.mt = b1*self.mt + (1-b1)*self.grad#*self.coef_grad
            self.vt = b2*self.vt + (1-b2)*np.sum(self.grad*self.grad)#*self.coef_grad**2
            mt_bias=1/(1-b1**(self.settings["step"]+1))*self.mt
            vt_bias=1/(1-b2**(self.settings["step"]+1))*self.vt
            
            if(self.settings["step"]>0):
                v_c_1=-eta*(vt_bias**(-0.5)+eps) * mt_bias
                vec_chang=v_c_1 +0.3*(v_c_1+ eta * (self.vt_bias_prev**(-0.5)+eps)*self.mt_bias_prev)
            else:
                vec_chang=-eta*(vt_bias+eps)**(-0.5) * mt_bias

            norm_chang=np.linalg.norm(vec_chang)
            if(norm_chang>TRUST_RAD):
                vec_chang=TRUST_RAD/norm_chang*vec_chang
            self.xyzs=self.xyzs+vec_chang
            
            self.vt_bias_prev=vt_bias
            self.mt_bias_prev=mt_bias
            
            self.update_xyzs_strs()
            self.Method.grad("!result")
            self.Method.read_grad() 
            
        elif 0:
            if(self.settings["step"]<1):
                #adam (rotational correction)
                self.prev_grad=copy.deepcopy(self.grad)
                
                self.apply_grad()
                self.update_xyzs_strs()
                self.Method.grad("!result")
                self.Method.read_grad()
            else:
                self.mt = b1*self.mt + (1-b1)*self.grad
                self.vt = b2*self.vt + (1-b2)*np.sum(self.grad*self.grad)

                cur_diff_mt = -self.prev_grad+self.grad
                self.diff_mt = b1_diff*self.diff_mt + (1-b1_diff)*cur_diff_mt
                self.diff_vt = b2*self.diff_vt + (1-b2)*np.sum(cur_diff_mt*cur_diff_mt)

                mt_bias=1/(1-b1**(self.settings["step"]))*self.mt
                vt_bias=1/(1-b2**(self.settings["step"]))*self.vt

                diff_mt_bias=1/(1-b1_diff**(self.settings["step"]))*self.diff_mt
                diff_vt_bias=1/(1-b2**(self.settings["step"]))*self.diff_vt

                vec_chang=-eta*(vt_bias**(-0.5)+eps) * mt_bias - eta_diff * mgsin_sqr * (diff_vt_bias**(-0.5)+eps) * diff_mt_bias 
                norm_chang=np.linalg.norm(vec_chang)
                if(norm_chang>TRUST_RAD):
                    vec_chang=TRUST_RAD/norm_chang*vec_chang
                
                self.xyzs=self.xyzs+vec_chang
                self.prev_grad=copy.deepcopy(self.grad)

                self.update_xyzs_strs()
                self.Method.grad("!result")
                self.Method.read_grad() 
        elif 0:
            init_xyzs=copy.deepcopy(self.xyzs)
            #optimistic ADAM
            self.mt = b1*self.mt + (1-b1)*self.grad
            self.vt = b2*self.vt + (1-b2)*np.sum(self.grad*self.grad)
            mt_bias=1/(1-b1**(self.settings["step"]+1))*self.mt
            vt_bias=1/(1-b2**(self.settings["step"]+1))*self.vt
            
            vec_chang=-eta*(vt_bias**(-0.5)+eps) * mt_bias
            norm_chang=np.linalg.norm(vec_chang)
            if(norm_chang>TRUST_RAD):
                vec_chang=TRUST_RAD/norm_chang*vec_chang
            self.xyzs=self.xyzs+vec_chang
            #rotating correction
            
            self.update_xyzs_strs()
            self.Method.grad("!result")
            self.Method.read_grad() 
            self.grad, maxgrad=self.get_grad()
            self.mirror()
            self.alter_grad()
            
            self.mt_r = b1*self.mt_r + (1-b1)*self.grad
            self.vt_r = b2*self.vt_r + (1-b2)*np.sum(self.grad*self.grad)
            
            mt_r_bias=1/(1-b1**(self.settings["step"]+1))*self.mt_r
            vt_r_bias=1/(1-b2**(self.settings["step"]+1))*self.vt_r

            vec_chang=-eta_r*(vt_r_bias**(-0.5)+eps) * mt_r_bias
            norm_chang=np.linalg.norm(vec_chang)
            if(norm_chang>TRUST_RAD):
                vec_chang=TRUST_RAD/norm_chang*vec_chang

            self.xyzs==init_xyzs+vec_chang

            self.update_xyzs_strs()
            self.Method.grad("!result")
            self.Method.read_grad() 
        else:#GD
            self.apply_grad()
            self.update_xyzs_strs()
            self.Method.grad("!result")
            self.Method.read_grad()
        
        self.settings["step"]+=1
        '''if maxgrad<self.prev_maxgrad:
            self.coef_grad*=1.01
        elif maxgrad>self.prev_maxgrad:
            self.coef_grad*=0.9
            if self.coef_grad<0.4:
                self.coef_grad=0.4
        
        print(f"coef grad {self.coef_grad}")'''
        self.prev_maxgrad=maxgrad
        return maxgrad
        
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
          mirror_coef=1,
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
    optTS(xyz_path=os.path.join("tests","piece_s_test", "to_opt.xyz"), threshold_rel=8, threshold_force=0.00001, mirror_coef=0.4, print_output=True, maxstep=10**4, programm=dict(name="xtb", force_constant= 6, acc=0.01),do_preopt=True,step_along=0)
    #optTS(xyz_path=os.path.join("tests","piece_s_test", "to_opt.xyz"), threshold_rel=8, threshold_force=0.00001, mirror_coef=0.4, print_output=True, maxstep=10**4, programm=dict(name="orca", memory="2000", nprocs=8, ORCA_PATH="/opt", method_str="BP86 def2-SVP"),do_preopt=True,step_along=0)
    '''
