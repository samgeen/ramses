'''
Make difference files for each code file in a folder
Pretty sure this is slower than doing it by hand but EH
Sam Geen, November 2013
'''

import os
import numpy as np

def run():
    mhds = list()
    files = os.listdir(".")
    for f1 in files:
        if ".f90" in f1:
            stub = f1[:-4]
            print stub
            f2 = "../../../mhd/"+f1
            if os.path.isfile(f2):
                mhds.append(f1)
            else:
                f2 = "../../../*/"+f1
            os.system("sdiff "+f1+" "+f2+" > diffs/diff_"+stub+".dat")
            np.savetxt("diffs/mhds.txt",mhds,delimiter="\n",fmt="%s")
            
       

if __name__=="__main__":
    run()
