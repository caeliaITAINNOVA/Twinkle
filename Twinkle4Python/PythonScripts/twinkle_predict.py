# -*- coding: utf-8 -*-

'''
TWINKLE python prediction

Version: 1.0
copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (Spain)
Date: 24 January 2020
Author: Instituto Tecnológico de Aragón (caelia@itainnova.es) - I.Viejo
All rights reserved	


TWINKLE LICENSE
Twinkle has two kind of licenses: commercial and open-source.

Commercial license
If you want to use Twinkle in a commercial way (any use out of open source applications under GNU GPL license v3 - https://www.gnu.org/licenses/gpl-3.0.html), please contact with ITAINNOVA (http://www.itainnova.es).

Open source license
If you are creating an open source application under a license compatible with the GNU GPL license v3, you may use Twinkle under the terms of the GPLv3 (https://www.gnu.org/licenses/gpl-3.0.html).



EXAMPLE OF USE:
filename='Results_file.txt'

values=[400,30,2]
prediction=twinkle_predict.predictwinkle(filename,values)


'''

#load global modules
import numpy as np

#load local modules


###################################
# Function readtwinkle
def predictwinkle(filename, params, **kwargs):
    '''
    Predict the value of params, if the value of the parameter is outside return the clossets one
        INPUTS:
            -filename (string): name of the file generated with Twinkle
            -params (list of list or numpy): rows the cases, columns are the value of each dimension of the case
        OUTPUTS:
            -predic (float or numpy): is the prediction
    '''
    
    with open(filename,'r') as f:
        lines=f.readlines()
    
    n_terms=int([line for i,line in enumerate(lines) if line.startswith('*** Number of terms:')][0].replace('\n','').split('\t')[-1])
    n_dims=int([line for i,line in enumerate(lines) if line.startswith('*** Number of dimensions:')][0].replace('\n','').split('\t')[-1])
    
    idx_discret=[i for i,line in enumerate(lines) if line.startswith('*** Dimension ')]
    
    discret=[]
    for n_dim in range(n_dims):
        tmp=lines[idx_discret[n_dim]].replace('\n','').split('\t')[1:n_terms+1]
        discret.append([int(v) for v in tmp])
    
    idx_terms=[i for i,line in enumerate(lines) if line.startswith('Term')]
    idx_alphas=[i for i,line in enumerate(lines) if line.startswith('Alpha')]
    idx_dimension=[i for i,line in enumerate(lines) if line.startswith('Dimension')]
 
    
    cases=np.asarray(params)
    if len(cases.shape)==1:
        cases=np.reshape(cases,(1,cases.shape[0]))
    
    # PREDICTION
    predic=0.0
    for n_term in range(n_terms):
        tmp_term=float(lines[idx_alphas[n_term]].replace('\n','').split('\t')[-1])
        for n_dim in range(n_dims):
            vec=lines[idx_dimension[n_term*n_dims+n_dim]+1:idx_dimension[n_term*n_dims+n_dim]+1+discret[n_dim][n_term]]
            values=np.asarray([v.replace('\n','').split('\t') for v in vec], dtype=float)
            tmp=np.interp(cases[:,n_dim],values[:,0],values[:,1])
            tmp_term=tmp_term*tmp
            
        predic+=tmp_term
    
    if len(predic)==1:
        return float(predic)
    else:
        return predic
