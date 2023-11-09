#! /usr/bin/python
### Protein check input function
#Packages
import numpy as np
import pandas as pd

def Protein_check(df,opt="other"):
    seqs=df['Sequence'].values
    ind=0
    std = list("ACDEFGHIKLMNPQRSTVWY")
    for i in range(len(seqs[0])):
        if not (seqs[0][i] in std):
            ind = i
            break
    if ind!=0:
        string='Error : Unexpected character in protein sequence in position: '+str(ind)+' '+str(seqs[0][ind])
        return string
    else:
        if opt=="B_cell":
            string='Correct format of the selected protein for the analysis. You can proceed either to Parameters Tab Panel to configure a custom analysis or you can proceed to  Execute panel from the menu to continue with the default options.'
        else:
            string='Correct format of the selected protein for the analysis. You can must proceed to Parameters Tab Panel in the menu to select the alleles for the epitope prediction.'
        return string

        
