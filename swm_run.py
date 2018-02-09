## SHALLOW WATER MODEL
from __future__ import print_function       # tested with python 3.6 and 2.7.12

from swm_param import set_param
from swm_integration import time_integration

## import all functions

# exec(open(path+'swm_param.py').read())
# exec(open(path+'swm_operators.py').read())
#exec(open(path+'swm_rhs.py').read())
#exec(open(path+'swm_integration.py').read())
#exec(open(path+'swm_output.py').read())

## set all parameters and initial conditions, and run model
u,v,eta = set_param()
u,v,eta = time_integration(u,v,eta)
