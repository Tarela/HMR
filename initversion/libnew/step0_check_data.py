#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *
import platform

# --------------------------
# custom package
# --------------------------


### tool function
from HMRpipe.Utility        import (sp,
                                   sperr,
                                   pdf_name,
                                   raise_error,
                                   detect_memory,
                                   wlog,
                                   ewlog,
                                   checkbedformat,
                                   CMD
                                   )
# --------------------------
# main 
# --------------------------

def step0_check_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    ### check data path , format ,
    if "~" in conf_dict['General']['HMRpeak']:
        ewlog('require absolute path for HMRpeak bed file, HMRpeak file cannot contain "~", current HMRpeak file is %s'%(conf_dict['General']['HMRpeak']),logfile)
    if "~" in conf_dict['General']['signal']:
        ewlog('require absolute path for HMsignal bigwig file, signal file cannot contain "~", current signal file is %s'%(conf_dict['General']['signal']),logfile)
    if not conf_dict['General']['HMRpeak'].startswith('/'):
        conf_dict['General']['HMRpeak'] = conf_dict['General']['startdir'] + conf_dict['General']['HMRpeak']
    if not conf_dict['General']['signal'].startswith('/'):
        conf_dict['General']['signal'] = conf_dict['General']['startdir'] + conf_dict['General']['signal']
    
    if not os.path.isfile(conf_dict['General']['HMRpeak']):
        ewlog("HMRpeak file %s not found"%(conf_dict['General']['HMRpeak']),logfile)
    if not os.path.isfile(conf_dict['General']['signal']):
        ewlog("signal bw file %s not found"%(conf_dict['General']['signal']),logfile)
        

    if not conf_dict['General']['HMRpeak'].endswith('.bed') :
        ewlog('extenion of HMR peak file is not .bed',logfile)
    checkbed = checkbedformat(conf_dict['General']['HMRpeak'],1000)
    if checkbed == "fail":
        ewlog("HMRpeak file is not a bed file",logfile)
    elif checkbed == "lesspeak":
        ewlog("HMRpeak file contains less than 1000 peaks") 

    if conf_dict['General']['signal'].endswith('.bw') :
        conf_dict['General']['signalname'] = conf_dict['General']['signal'].split("/")[-1][:-3]
    elif conf_dict['General']['signal'].endswith('.bigwig'):
        conf_dict['General']['signalname'] = conf_dict['General']['signal'].split("/")[-1][:-7]
    else:
        wlog('[WARNING] extension of signal bw file is not bw/bigwig',logfile)

    ### check TFpeak folder    
    if "~" in conf_dict['General']['peakFolder']:
        ewlog('require absolute path for peakFolder, peakFolder cannot contain "~", current peakFolder is %s'%(conf_dict['General']['peakFolder']),logfile)
    if not conf_dict['General']['peakFolder'].startswith('/'):
        conf_dict['General']['peakFolder'] = conf_dict['General']['startdir'] + conf_dict['General']['peakFolder']
    if not conf_dict['General']['peakFolder'].endswith('/'):
        conf_dict['General']['peakFolder'] += "/"
    if not os.path.isdir(conf_dict['General']['peakFolder']):
        ewlog("peakFolder %s not found"%(conf_dict['General']['peakFolder']),logfile)

    wlog("Check the peak.bed files in the peakFolder, only '.bed' files with >1000 peaks are included in the following analysis",logfile)
    conf_dict['General']['peakfilenames'] = []
    for f in os.listdir(conf_dict['General']['peakFolder']):
        if f.endswith(".bed") and os.path.isfile(conf_dict['General']['peakFolder']+f):
            checkbed = checkbedformat(conf_dict['General']['HMRpeak'],1000)
            if checkbed == "pass":
                conf_dict['General']['peakfilenames'].append(f[:-4])
    if (len(conf_dict['General']['peakfilenames']) == 0):
        ewlog("no peak file (cofactor candidate) in bed format & >1000peaks are included, exit",logfile)
    else:
        wlog("%s peak files (cofactor candidates) are included"%(len(conf_dict['General']['peakfilenames'])),logfile)

    outf = open(conf_dict['General']['outname']+"_cofactor_candidate_list.txt",'w')
    for cofactor in conf_dict['General']['peakfilenames']:
        outf.write(cofactor+"\n")
    outf.close()
    ### check options
    wlog('check option: ',logfile)

    try:
        wlog("extend length is %s bp"%(int(conf_dict['options']['ext'])),logfile)
        conf_dict['options']['ext']=int(conf_dict['options']['ext'])
    except:
        wlog("extend length %s is not valid, use default value: 1000bp"%(conf_dict['options']['ext']),logfile)
        conf_dict['options']['ext'] = 1000

    if conf_dict['options']['Qvalue']:
        conf_dict['options']['usePQ'] = "Q"
        try:
            wlog("use Qvalue = %s as cutoff"%(str(float(conf_dict['options']['Qvalue']))),logfile)
        except:
            wlog("input Qvalue %s is not recognized, use default Qvalue=0.01"%(conf_dict['options']['Qvalue']),logfile)
            conf_dict['options']['Qvalue'] = 0.01 
        if float(conf_dict['options']['Qvalue']) >= 1:
            wlog("input Qvalue %s is not valid, use default Qvalue=0.01"%(conf_dict['options']['Qvalue']),logfile)
            conf_dict['options']['Qvalue'] = 0.01 
    elif conf_dict['options']['Pvalue']:
        conf_dict['options']['usePQ'] = "P"
        try:
            wlog("use Pvalue = %s as cutoff"%(str(float(conf_dict['options']['Pvalue']))),logfile)
        except:
            wlog("input Pvalue %s is not recognized, use default Pvalue=0.001"%(conf_dict['options']['Pvalue']),logfile)
            conf_dict['options']['Pvalue'] = 0.001 
        if float(conf_dict['options']['Pvalue']) >= 1:
            wlog("input Pvalue %s is not valid, use default Pvalue=0.001"%(conf_dict['options']['Pvalue']),logfile)
            conf_dict['options']['Pvalue'] = 0.001 
    else:
        wlog("use Qvalue = 0.01 as cutoff",logfile)
        conf_dict['options']['usePQ'] = "Q"
        conf_dict['options']['Qvalue'] = 0.01

    if conf_dict['options']['usePQ']:
        conf_dict['options']['cutoff']=conf_dict['options']['Qvalue']
    else:
        conf_dict['options']['cutoff']=conf_dict['options']['Pvalue']

    try: 
        usealpha = float(conf_dict['options']['Alpha'])
        if usealpha >=1:
            wlog("alpha cannot be >=1, use alpha=0.5",logfile)
            conf_dict['options']['Alpha'] = 0.5
        else:
            wlog("Alpha = %s"(str(float(conf_dict['options']['Alpha']))),logfile)
            conf_dict['options']['Alpha']=usealpha
    except:
        wlog("input alpha %s is not valid, use alpha=0.5"%(conf_dict['options']['Alpha']),logfile)
    
    wlog("Lambda choice is %s"%(conf_dict['options']['Lambda']),logfile)
    if conf_dict['options']['TopNcofactors'] == "all":
        wlog("all significant co-factors will be output",logfile)
    else:
        try:
            topTF = int(conf_dict['options']['TopNcofactors'])
            wlog("the topN number %s will be output"%(conf_dict['options']['TopNcofactors']),logfile)
            conf_dict['options']['TopNcofactors'] = topTF
        except:
            wlog("the topN number %s is not valid, output top5 co-factors"%(conf_dict['options']['TopNcofactors']),logfile)
            conf_dict['options']['TopNcofactors'] = 5

    OS = platform.system()
    if OS == "Linux":
        bwsum_software = "bigWigSummary_linux"
    elif OS == "Darwin":
        bwsum_software = "bigWigSummary_mac"
    else:
        wlog("detected system is nither linux nor mac, try linux version",logfile)
        bwsum_software = "bigWigSummary_linux"

    conf_dict['General']['software'] = bwsum_software
    ### check Rscript
    #if not 'Usage' in sperr('Rscript')[1] and not 'version' in sperr('Rscript')[1]:
    #    ewlog('require Rscript',logfile)
    
    ### check pdflatex
    if sp('pdflatex --help')[0] == "":
        wlog('pdflatex was not installed, HMR is still processing but no summary report generated',logfile)
        conf_dict['General']['latex'] = 0
    else:
        conf_dict['General']['latex'] = 1

    return conf_dict
    
    
    
