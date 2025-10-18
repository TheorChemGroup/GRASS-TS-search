import matplotlib.pyplot as plt
import numpy as np

colors_through=("blue","red","green","orange")
species = (0,0.2,	0.4,	0.6,	0.8,	1,	1.5,	2)
values = {
    "Mirror": (82,	5,	5,	10,	4,	12,	2,	1),
    'Neb-TS': (65,	6,	8,	7,	7,	13,	11,	2),
    'NT2': (32,	6,	10,	23,	18,	22,	5,	2),
    'preopt': (29,	0,	6,	10,	11,	41,	21,	3),
}

x = np.arange(len(species))  # the label locations
width = 1/(len(values)+1)  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained',figsize=(10, 6))

for attribute, measurement in values.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    labels=ax.bar_label(rects, padding=2,rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('Testing of existing methods on GFN1 level of theory on zba121 testset',size=16)

ax.set_xlabel("RMSD",size=14)
ax.set_xticks(x - width, species,size=12)

ax.set_ylabel('Number of optimized with RMSD within interval',size=14)
ax.set_ylim(0, 92)
ax.set_yticks(range(0,91,10),range(0,91,10),size=12)

ax.legend(["Mirror", "Neb-TS", "NT2","preopt"],loc=(0.885,0.35))




ax_corner = fig.add_axes([0.5, 0.6, 0.4955, 0.3])

import numpy as np
import pandas as pd

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

    sig_sq=70
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


ax_corner.set_title('Probability density function')

x=np.linspace(0,2.5,500)
for name in ("Mirror","Neb-TS","NT2","preopt"):
    __f=lambda x_i:den_f(x_i,name)
    print(integrate(__f,0,3,20))
    ax_corner.plot(x,np.vectorize(den_f)(x,name),label=name,linewidth=2)



ax_corner.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True) # Hide ticks and labels
ax_corner.set_xlabel("RMSD",labelpad=-5)
#ax_corner.legend()

#fig.legend(loc='upper right', ncols=4)


fig.savefig(f"pictures/graph", dpi=300)




plt.savefig(f"pictures/chart", dpi=300)