from TS_find import optTS
import sys
import os

path=""
mode="autostrict"
ratio=0
optimized_cap=0
print_output=False
if (True if len(sys.argv)<2 else sys.argv[1]=="-h"):
    print(f"usage: {sys.argv[0]} [--path <path>] [-h] [--ratio <ratio>] [--opt_cap <opt_cap>] [-p] [--mode <mode>] [--steps <max steps>]")
    exit()
xyz_name=sys.argv[1]

i=2
while i < len(sys.argv):
    if sys.argv[i]=="--path":
        path=sys.argv[i+1]
        i+=1
    elif sys.argv[i]=="--ratio":
        ratio=float(sys.argv[i+1])
        i+=1
    elif sys.argv[i]=="--opt_cap":
        optimized_cap=float(sys.argv[i+1])
        i+=1
    elif sys.argv[i]=="-p":
        print_output=True
    elif sys.argv[i]=="--mode":
        mode=sys.argv[i+1]
        i+=1
    elif sys.argv[i]=="--steps":
        maxstep=int(sys.argv[i+1])
        i+=1
    else:
        print(f"unknown parameter {sys.argv[i]} at position {i}")
        exit()
    i+=1

optTS(os.path.join(os.getcwd(),path),xyz_name=xyz_name, ratio=ratio, optimized_cap=optimized_cap, print_output=True,mode=mode)