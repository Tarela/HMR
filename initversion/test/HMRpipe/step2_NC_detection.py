#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time

# --------------------------
# custom package
# --------------------------

### tool function
from HMRpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   createDIR)
# --------------------------
# main 
# --------------------------
def step2_NC_detection(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''   
    # start
    # create section for 
    # Rscript detectNonCanonical.r outname signalname usePQ cutoff alpha lambdachoice topN
    cmd = "Rscript %s %s %s %s %s %s %s %s"%(conf_dict['rscript']+"detectNonCanonical.r",
                         conf_dict['General']['outname'],
                         conf_dict['General']['signalname'],
                         conf_dict['options']['usePQ'],
                         conf_dict['options']['cutoff'],
                         conf_dict['options']['Alpha'],
                         conf_dict['options']['Lambda'],
                         conf_dict['options']['TopNcofactors'])
    rwlog(cmd,logfile)

    return conf_dict


