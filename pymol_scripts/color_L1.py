#
# -- basicCodeBlock.py
#
import numpy as np
from pymol import cmd, stored

def color_L1( arg1, arg2, arg3, arg4=1 ):
    '''
DESCRIPTION

    Brief description: Colors the identified L1 features 
    '''
    # Running on pymol cmdline
    #run ./pymol_scripts/color_L1.py
    #color_L1 78 137 158 188 241 292 573, green, 2
    
    Feature_set=[int(j) for j in np.array(arg1.split())]    
    col=arg2
    L=int(arg3)

    #reset color to gray
    if arg4==1:
        cmd.color('gray')
    
    F=np.array([np.arange(i-L,i+L+1) for i in Feature_set]).flatten()    
    for feature in F:
        print("resi "+str(feature))
        cmd.color(col, "resi "+ str(feature))
    

cmd.extend( "color_L1", color_L1 );
