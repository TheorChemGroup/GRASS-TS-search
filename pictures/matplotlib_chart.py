import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

test_f = lambda x:np.sin(x)

def integrate(f,a,b,n):
    res=0
    
    lb=a
    flb=f(lb)

    l_ab=b-a
    rev_n=1/n
    ab_rev_n_d2=l_ab*0.5/n
    for i in range(1,n+1):
        rb=a+l_ab*i*rev_n
        frb=f(rb)
        res+=(flb+frb)*ab_rev_n_d2
        flb=frb
        lb=rb

    return res
        

def den_f(x,name):
    res=0

    sig_sq=30
    norm_coef=(sig_sq/(2*np.pi))**0.5
    
    cnt=0
    for val in df[name]:
        if(val==val):
            cnt+=1
    
    for val in df[name]:
        if(val==val):
            res+=0.5*norm_coef*np.exp(-0.5*sig_sq*(val-x)**2)
            res+=0.5*norm_coef*np.exp(-0.5*sig_sq*(val+x)**2)#продление в отрицательную часть, чтобы в нуле было ровно
    return 2*res/cnt
            
df=pd.read_csv("pictures/rmsds_new2.csv")

fig, ax_main = plt.subplots(figsize=(10, 6))

x=np.linspace(0,3,500)
for name in ("result", "preopt", "nt2","neb2"):
    __f=lambda x_i:den_f(x_i,name)
    print(integrate(__f,0,10,100))
    ax_main.plot(x,np.vectorize(den_f)(x,name),label=name)

ax_corner = fig.add_axes([0.6, 0.6, 0.2, 0.2])
x_corner = np.linspace(0, 5, 50)
y_corner = np.cos(x_corner)
ax_corner.plot(x_corner, y_corner, 'r--', label='Corner Plot Data')
ax_corner.set_title('Inset Plot')
ax_corner.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False) # Hide ticks and labels
ax_corner.legend()

fig.legend(loc='upper right', ncols=4)


fig.savefig(f"pictures/graph", dpi=300)